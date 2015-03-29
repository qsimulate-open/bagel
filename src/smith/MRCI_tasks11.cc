//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks11.cc
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

#include <src/smith/MRCI_tasks11.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task500::Task_local::compute() {
  const Index c2 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index x2 = b(3);
  // tensor label: I82
  std::unique_ptr<double[]> odata = out()->move_block(c2, x3, a1, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x3, a1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x3, a1, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c3, a4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c3.size(), a4.size());
      // tensor label: I823
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x2, c2, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x2, c2, c3)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a4.size(), x2.size(), c2.size(), c3.size());
      dgemm_("T", "N", x3.size()*a1.size(), x2.size()*c2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, x3.size()*a1.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size(), x2.size(), c2.size());
  out()->put_block(odata, c2, x3, a1, x2);
}

void Task501::Task_local::compute() {
  const Index a4 = b(0);
  const Index x2 = b(1);
  const Index c2 = b(2);
  const Index c3 = b(3);
  // tensor label: I823
  std::unique_ptr<double[]> odata = out()->move_block(a4, x2, c2, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x2, c2, c3);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, a4.size(), x2.size(), c2.size(), c3.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c2, x2, a4, c3);
    sort_indices<2,1,0,3,1,1,2,1>(i1data, odata, c2.size(), x2.size(), a4.size(), c3.size());
  }
  out()->put_block(odata, a4, x2, c2, c3);
}

void Task502::Task_local::compute() {
  const Index c2 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index x2 = b(3);
  // tensor label: I82
  std::unique_ptr<double[]> odata = out()->move_block(c2, x3, a1, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x3, a1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x3, a1, x2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c3, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a4, c3, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c3.size(), a1.size());
      // tensor label: I826
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x2, a4, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x2, a4, c3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), x2.size(), a4.size(), c3.size());
      dgemm_("T", "N", x3.size()*a1.size(), c2.size()*x2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, x3.size()*a1.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size(), c2.size(), x2.size());
  out()->put_block(odata, c2, x3, a1, x2);
}

void Task503::Task_local::compute() {
  const Index c2 = b(0);
  const Index x2 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I826
  std::unique_ptr<double[]> odata = out()->move_block(c2, x2, a4, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x2, a4, c3);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c2.size(), x2.size(), a4.size(), c3.size());
  }
  out()->put_block(odata, c2, x2, a4, c3);
}

void Task504::Task_local::compute() {
  const Index c2 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index x2 = b(3);
  // tensor label: I82
  std::unique_ptr<double[]> odata = out()->move_block(c2, x3, a1, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x3, a1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x3, a1, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a4.size());
      // tensor label: I835
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x2, a3, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x2, a3, a1)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), x2.size(), a3.size(), a1.size());
      dgemm_("T", "N", x3.size()*c2.size(), x2.size()*a1.size(), a4.size()*a3.size(),
             1.0, i0data_sorted, a4.size()*a3.size(), i1data_sorted, a4.size()*a3.size(),
             1.0, odata_sorted, x3.size()*c2.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x3.size(), c2.size(), x2.size(), a1.size());
  out()->put_block(odata, c2, x3, a1, x2);
}

void Task505::Task_local::compute() {
  const Index a4 = b(0);
  const Index x2 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I835
  std::unique_ptr<double[]> odata = out()->move_block(a4, x2, a3, a1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x2, a3, a1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, a4.size(), x2.size(), a3.size(), a1.size());
  }
  out()->put_block(odata, a4, x2, a3, a1);
}

void Task506::Task_local::compute() {
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
      // tensor label: Gamma29
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I88
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x3.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), c2.size()*a1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a1.size());
  out()->put_block(odata, c2, x1, x0, a1);
}

void Task507::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I88
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x3, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a1, x3, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a1, x3, x2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a1.size(), x3.size(), x2.size());
    // tensor label: I89
    std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c2.size(), c3.size());
    dgemm_("T", "N", a1.size()*x3.size()*x2.size(), c2.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, a1.size()*x3.size()*x2.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a1.size(), x3.size(), x2.size(), c2.size());
  out()->put_block(odata, c2, a1, x3, x2);
}

void Task508::Task_local::compute() {
  const Index c2 = b(0);
  const Index c3 = b(1);
  // tensor label: I89
  std::unique_ptr<double[]> odata = out()->move_block(c2, c3);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c3);
    sort_indices<0,1,1,1,1,1>(i0data, odata, c2.size(), c3.size());
  }
  out()->put_block(odata, c2, c3);
}

void Task509::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I88
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x3, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, x3, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, x3, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), x3.size(), x2.size());
    // tensor label: I92
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, a1)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), a1.size());
    dgemm_("T", "N", c2.size()*x3.size()*x2.size(), a1.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, c2.size()*x3.size()*x2.size());
  }
  sort_indices<0,3,1,2,1,1,1,1>(odata_sorted, odata, c2.size(), x3.size(), x2.size(), a1.size());
  out()->put_block(odata, c2, a1, x3, x2);
}

void Task510::Task_local::compute() {
  const Index a3 = b(0);
  const Index a1 = b(1);
  // tensor label: I92
  std::unique_ptr<double[]> odata = out()->move_block(a3, a1);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a1);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, a3.size(), a1.size());
  }
  out()->put_block(odata, a3, a1);
}

void Task511::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I88
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x3, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a1.size());
    // tensor label: I104
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x2)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), x2.size());
    dgemm_("T", "N", x3.size()*c2.size()*a1.size(), x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x3.size()*c2.size()*a1.size());
  }
  sort_indices<1,2,0,3,1,1,1,1>(odata_sorted, odata, x3.size(), c2.size(), a1.size(), x2.size());
  out()->put_block(odata, c2, a1, x3, x2);
}

void Task512::Task_local::compute() {
  const Index a3 = b(0);
  const Index x2 = b(1);
  // tensor label: I104
  std::unique_ptr<double[]> odata = out()->move_block(a3, x2);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x2);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, a3.size(), x2.size());
  }
  out()->put_block(odata, a3, x2);
}

void Task513::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I88
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x3, x2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a4, c3, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c3.size(), x2.size());
      // tensor label: I727
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, a4, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, a4, c3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), a4.size(), c3.size());
      dgemm_("T", "N", x3.size()*x2.size(), c2.size()*a1.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c2.size(), a1.size());
  out()->put_block(odata, c2, a1, x3, x2);
}

void Task514::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I727
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, a4, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, a4, c3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c2.size(), a1.size(), a4.size(), c3.size());
  }
  out()->put_block(odata, c2, a1, a4, c3);
}

void Task515::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I88
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x3, x2), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a3, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a3, x3, x2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a3.size(), x3.size(), x2.size());
      // tensor label: I754
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c4, a3, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c4, a3, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), c4.size(), a3.size(), a1.size());
      dgemm_("T", "N", x3.size()*x2.size(), c2.size()*a1.size(), c4.size()*a3.size(),
             1.0, i0data_sorted, c4.size()*a3.size(), i1data_sorted, c4.size()*a3.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c2.size(), a1.size());
  out()->put_block(odata, c2, a1, x3, x2);
}

void Task516::Task_local::compute() {
  const Index c2 = b(0);
  const Index c4 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I754
  std::unique_ptr<double[]> odata = out()->move_block(c2, c4, a3, a1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c4, a3, a1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c2.size(), c4.size(), a3.size(), a1.size());
  }
  out()->put_block(odata, c2, c4, a3, a1);
}

void Task517::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I88
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x3, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4, x3, x2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), x3.size(), x2.size());
      // tensor label: I757
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, a4, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, a4, c3)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), a4.size(), c3.size());
      dgemm_("T", "N", x3.size()*x2.size(), c2.size()*a1.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c2.size(), a1.size());
  out()->put_block(odata, c2, a1, x3, x2);
}

void Task518::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I757
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, a4, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, a4, c3);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c2.size(), a1.size(), a4.size(), c3.size());
  }
  out()->put_block(odata, c2, a1, a4, c3);
}

void Task519::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I88
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x3, x2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a4, c3, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a4, c3, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a4.size(), c3.size(), a1.size());
      // tensor label: I772
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3, x3, x2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size(), x3.size(), x2.size());
      dgemm_("T", "N", c2.size()*a1.size(), x3.size()*x2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x3.size(), x2.size());
  out()->put_block(odata, c2, a1, x3, x2);
}

void Task520::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I772
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3, x3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3, x3, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, a4.size(), c3.size(), x3.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x3, x2, a4, c3);
    sort_indices<2,3,0,1,1,1,1,1>(i1data, odata, x3.size(), x2.size(), a4.size(), c3.size());
  }
  out()->put_block(odata, a4, c3, x3, x2);
}

void Task521::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I88
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x3, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c3, a4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c3.size(), a4.size());
      // tensor label: I775
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3, x3, x2)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size(), x3.size(), x2.size());
      dgemm_("T", "N", c2.size()*a1.size(), x3.size()*x2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x3.size(), x2.size());
  out()->put_block(odata, c2, a1, x3, x2);
}

void Task522::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I775
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

void Task523::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I88
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x3, x2), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c4, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c4, a3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c4.size(), a3.size());
      // tensor label: I787
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c4, a3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c4, a3, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), c4.size(), a3.size(), x2.size());
      dgemm_("T", "N", c2.size()*a1.size(), x3.size()*x2.size(), c4.size()*a3.size(),
             1.0, i0data_sorted, c4.size()*a3.size(), i1data_sorted, c4.size()*a3.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x3.size(), x2.size());
  out()->put_block(odata, c2, a1, x3, x2);
}

void Task524::Task_local::compute() {
  const Index x3 = b(0);
  const Index c4 = b(1);
  const Index a3 = b(2);
  const Index x2 = b(3);
  // tensor label: I787
  std::unique_ptr<double[]> odata = out()->move_block(x3, c4, a3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, c4, a3, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), c4.size(), a3.size(), x2.size());
  }
  out()->put_block(odata, x3, c4, a3, x2);
}

void Task525::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I88
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x3, x2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c3, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a4, c3, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c3.size(), a1.size());
      // tensor label: I820
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x2, c2, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x2, c2, c3)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, a4.size(), x2.size(), c2.size(), c3.size());
      dgemm_("T", "N", x3.size()*a1.size(), x2.size()*c2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, x3.size()*a1.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size(), x2.size(), c2.size());
  out()->put_block(odata, c2, a1, x3, x2);
}

void Task526::Task_local::compute() {
  const Index a4 = b(0);
  const Index x2 = b(1);
  const Index c2 = b(2);
  const Index c3 = b(3);
  // tensor label: I820
  std::unique_ptr<double[]> odata = out()->move_block(a4, x2, c2, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x2, c2, c3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, a4.size(), x2.size(), c2.size(), c3.size());
  }
  out()->put_block(odata, a4, x2, c2, c3);
}

void Task527::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I88
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x3, x2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a4, c2, a3)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c2.size(), a3.size());
      // tensor label: I832
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x2, a3, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x2, a3, a1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), x2.size(), a3.size(), a1.size());
      dgemm_("T", "N", x3.size()*c2.size(), x2.size()*a1.size(), a4.size()*a3.size(),
             1.0, i0data_sorted, a4.size()*a3.size(), i1data_sorted, a4.size()*a3.size(),
             1.0, odata_sorted, x3.size()*c2.size());
    }
  }
  sort_indices<1,3,0,2,1,1,1,1>(odata_sorted, odata, x3.size(), c2.size(), x2.size(), a1.size());
  out()->put_block(odata, c2, a1, x3, x2);
}

void Task528::Task_local::compute() {
  const Index a4 = b(0);
  const Index x2 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I832
  std::unique_ptr<double[]> odata = out()->move_block(a4, x2, a3, a1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x2, a3, a1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, a4.size(), x2.size(), a3.size(), a1.size());
  }
  out()->put_block(odata, a4, x2, a3, a1);
}

void Task529::Task_local::compute() {
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
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), x2.size());
    // tensor label: I94
    std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x0, x1, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x0, x1, x2)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a1.size(), x0.size(), x1.size(), x2.size());
    dgemm_("T", "N", c2.size(), a1.size()*x0.size()*x1.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, c2.size());
  }
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x0.size(), x1.size());
  out()->put_block(odata, c2, x1, x0, a1);
}

void Task530::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: I94
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x1, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, x2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma31
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x1, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x1, x2)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x1.size(), x2.size());
        // tensor label: I95
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a1, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a1, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
        dgemm_("T", "N", x0.size()*x1.size()*x2.size(), a1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x1.size()*x2.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), a1.size());
  out()->put_block(odata, a1, x0, x1, x2);
}

void Task531::Task_local::compute() {
  const Index x5 = b(0);
  const Index a1 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I95
  std::unique_ptr<double[]> odata = out()->move_block(x5, a1, x4, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x5.size(), a1.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, x5, a1, x4, x3);
}

void Task532::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  // tensor label: Gamma32
  std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
  // tensor label: I97
  std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size());
  dgemm_("T", "N", x1.size()*x0.size(), c2.size()*a1.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, x1.size()*x0.size());
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a1.size());
  out()->put_block(odata, c2, x1, x0, a1);
}

void Task533::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  // tensor label: I97
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a4, c3, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a4, c3, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a4.size(), c3.size(), a1.size());
      // tensor label: I98
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size());
      dgemm_("T", "N", c2.size()*a1.size(), 1, a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size());
  out()->put_block(odata, c2, a1);
}

void Task534::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  // tensor label: I98
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3);
    sort_indices<0,1,1,1,2,1>(i0data, odata, a4.size(), c3.size());
  }
  out()->put_block(odata, a4, c3);
}

void Task535::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  // tensor label: I97
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c3, a4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c3.size(), a4.size());
      // tensor label: I101
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size());
      dgemm_("T", "N", c2.size()*a1.size(), 1, a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size());
  out()->put_block(odata, c2, a1);
}

void Task536::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  // tensor label: I101
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3);
    sort_indices<0,1,1,1,-4,1>(i0data, odata, a4.size(), c3.size());
  }
  out()->put_block(odata, a4, c3);
}

void Task537::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  // tensor label: I97
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      for (auto& c5 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, c5, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4, c5, a1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), c5.size(), a1.size());
        // tensor label: I796
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c5, a4, c3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c5, a4, c3)]);
        sort_indices<3,2,1,0,0,1,1,1>(i1data, i1data_sorted, c2.size(), c5.size(), a4.size(), c3.size());
        dgemm_("T", "N", a1.size(), c2.size(), c5.size()*a4.size()*c3.size(),
               1.0, i0data_sorted, c5.size()*a4.size()*c3.size(), i1data_sorted, c5.size()*a4.size()*c3.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size());
  out()->put_block(odata, c2, a1);
}

void Task538::Task_local::compute() {
  const Index c2 = b(0);
  const Index c5 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I796
  std::unique_ptr<double[]> odata = out()->move_block(c2, c5, a4, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c5, a4, c3);
    sort_indices<0,1,2,3,1,1,4,1>(i0data, odata, c2.size(), c5.size(), a4.size(), c3.size());
  }
  out()->put_block(odata, c2, c5, a4, c3);
}

void Task539::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  // tensor label: I97
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& c5 : *range_[0]) {
      for (auto& a4 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a1, c5, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a1, c5, a4)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a1.size(), c5.size(), a4.size());
        // tensor label: I799
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c5, a4, c3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c5, a4, c3)]);
        sort_indices<3,1,2,0,0,1,1,1>(i1data, i1data_sorted, c2.size(), c5.size(), a4.size(), c3.size());
        dgemm_("T", "N", a1.size(), c2.size(), c5.size()*a4.size()*c3.size(),
               1.0, i0data_sorted, c5.size()*a4.size()*c3.size(), i1data_sorted, c5.size()*a4.size()*c3.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size());
  out()->put_block(odata, c2, a1);
}

void Task540::Task_local::compute() {
  const Index c2 = b(0);
  const Index c5 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I799
  std::unique_ptr<double[]> odata = out()->move_block(c2, c5, a4, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c5, a4, c3);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c2.size(), c5.size(), a4.size(), c3.size());
  }
  out()->put_block(odata, c2, c5, a4, c3);
}

void Task541::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  // tensor label: I97
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1), 0.0);
  for (auto& a5 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      for (auto& a4 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a5, c3, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a5, c3, a4)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), a5.size(), c3.size(), a4.size());
        // tensor label: I802
        std::unique_ptr<double[]> i1data = in(1)->get_block(a5, a1, a4, c3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, a1, a4, c3)]);
        sort_indices<0,3,2,1,0,1,1,1>(i1data, i1data_sorted, a5.size(), a1.size(), a4.size(), c3.size());
        dgemm_("T", "N", c2.size(), a1.size(), a5.size()*a4.size()*c3.size(),
               1.0, i0data_sorted, a5.size()*a4.size()*c3.size(), i1data_sorted, a5.size()*a4.size()*c3.size(),
               1.0, odata_sorted, c2.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size());
  out()->put_block(odata, c2, a1);
}

void Task542::Task_local::compute() {
  const Index a5 = b(0);
  const Index a1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I802
  std::unique_ptr<double[]> odata = out()->move_block(a5, a1, a4, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, a1, a4, c3);
    sort_indices<0,1,2,3,1,1,-4,1>(i0data, odata, a5.size(), a1.size(), a4.size(), c3.size());
  }
  out()->put_block(odata, a5, a1, a4, c3);
}

void Task543::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  // tensor label: I97
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      for (auto& a5 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a4, c3, a5);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a4, c3, a5)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), a4.size(), c3.size(), a5.size());
        // tensor label: I805
        std::unique_ptr<double[]> i1data = in(1)->get_block(a5, a1, a4, c3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, a1, a4, c3)]);
        sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a5.size(), a1.size(), a4.size(), c3.size());
        dgemm_("T", "N", c2.size(), a1.size(), a5.size()*a4.size()*c3.size(),
               1.0, i0data_sorted, a5.size()*a4.size()*c3.size(), i1data_sorted, a5.size()*a4.size()*c3.size(),
               1.0, odata_sorted, c2.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size());
  out()->put_block(odata, c2, a1);
}

void Task544::Task_local::compute() {
  const Index a5 = b(0);
  const Index a1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I805
  std::unique_ptr<double[]> odata = out()->move_block(a5, a1, a4, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, a1, a4, c3);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, a5.size(), a1.size(), a4.size(), c3.size());
  }
  out()->put_block(odata, a5, a1, a4, c3);
}

void Task545::Task_local::compute() {
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
      for (auto& c3 : *range_[0]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, c3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, c3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), c3.size());
        // tensor label: I654
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c3, x1, x2, x3, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c3, x1, x2, x3, x0)]);
        sort_indices<4,3,1,0,2,5,0,1,1,1>(i1data, i1data_sorted, c2.size(), c3.size(), x1.size(), x2.size(), x3.size(), x0.size());
        dgemm_("T", "N", a1.size(), c2.size()*x1.size()*x0.size(), c3.size()*x2.size()*x3.size(),
               1.0, i0data_sorted, c3.size()*x2.size()*x3.size(), i1data_sorted, c3.size()*x2.size()*x3.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), x1.size(), x0.size());
  out()->put_block(odata, c2, x1, x0, a1);
}

void Task546::Task_local::compute() {
  const Index c2 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  const Index x3 = b(4);
  const Index x0 = b(5);
  // tensor label: I654
  std::unique_ptr<double[]> odata = out()->move_block(c2, c3, x1, x2, x3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c3, x1, x2, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c3, x1, x2, x3, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: Gamma215
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x5, x2, x4, x3, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x5, x2, x4, x3, x0)]);
      sort_indices<1,3,0,2,4,5,0,1,1,1>(i0data, i0data_sorted, x1.size(), x5.size(), x2.size(), x4.size(), x3.size(), x0.size());
      // tensor label: I655
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x5, c3, x4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x5, c3, x4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, c2.size(), x5.size(), c3.size(), x4.size());
      dgemm_("T", "N", x1.size()*x2.size()*x3.size()*x0.size(), c2.size()*c3.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, x1.size()*x2.size()*x3.size()*x0.size());
    }
  }
  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, x1.size(), x2.size(), x3.size(), x0.size(), c2.size(), c3.size());
  out()->put_block(odata, c2, c3, x1, x2, x3, x0);
}

void Task547::Task_local::compute() {
  const Index c2 = b(0);
  const Index x5 = b(1);
  const Index c3 = b(2);
  const Index x4 = b(3);
  // tensor label: I655
  std::unique_ptr<double[]> odata = out()->move_block(c2, x5, c3, x4);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x5, c3, x4);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, c2.size(), x5.size(), c3.size(), x4.size());
  }
  out()->put_block(odata, c2, x5, c3, x4);
}

void Task548::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x4, a1, x3, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x4, a1, x3, x2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x4.size(), a1.size(), x3.size(), x2.size());
        // tensor label: I657
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x1, x4, x0, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x1, x4, x0, x3, x2)]);
        sort_indices<2,4,5,0,1,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), x1.size(), x4.size(), x0.size(), x3.size(), x2.size());
        dgemm_("T", "N", a1.size(), c2.size()*x1.size()*x0.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), x1.size(), x0.size());
  out()->put_block(odata, c2, x1, x0, a1);
}

void Task549::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  const Index x3 = b(4);
  const Index x2 = b(5);
  // tensor label: I657
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, x4, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x4, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x4, x0, x3, x2), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: Gamma216
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x1, x5, x4, x0, x3, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, x1, x5, x4, x0, x3, x2)]);
        sort_indices<0,1,3,2,4,5,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x1.size(), x5.size(), x4.size(), x0.size(), x3.size(), x2.size());
        // tensor label: I658
        std::unique_ptr<double[]> i1data = in(1)->get_block(x7, x6, c2, x5);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x7, x6, c2, x5)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x7.size(), x6.size(), c2.size(), x5.size());
        dgemm_("T", "N", x1.size()*x4.size()*x0.size()*x3.size()*x2.size(), c2.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, x1.size()*x4.size()*x0.size()*x3.size()*x2.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x1.size(), x4.size(), x0.size(), x3.size(), x2.size(), c2.size());
  out()->put_block(odata, c2, x1, x4, x0, x3, x2);
}

#endif
