//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks6.cc
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

#include <src/smith/MRCI_tasks6.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task250::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2, c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c4, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c4, x1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c4.size(), x1.size());
    // tensor label: I41
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c4)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), c4.size());
    dgemm_("T", "N", c3.size()*a2.size()*x1.size(), c1.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c3.size()*a2.size()*x1.size());
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), x1.size(), c1.size());
  out()->put_block(odata, c1, a2, c3, x1);
}

void Task251::Task_local::compute() {
  const Index c1 = b(0);
  const Index c4 = b(1);
  // tensor label: I41
  std::unique_ptr<double[]> odata = out()->move_block(c1, c4);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c4);
    sort_indices<0,1,1,1,1,1>(i0data, odata, c1.size(), c4.size());
  }
  out()->put_block(odata, c1, c4);
}

void Task252::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2, c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c1, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a2, c1, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c1.size(), x1.size());
    // tensor label: I44
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, c4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, c4)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), c4.size());
    dgemm_("T", "N", a2.size()*c1.size()*x1.size(), c3.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, a2.size()*c1.size()*x1.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x1.size(), c3.size());
  out()->put_block(odata, c1, a2, c3, x1);
}

void Task253::Task_local::compute() {
  const Index c3 = b(0);
  const Index c4 = b(1);
  // tensor label: I44
  std::unique_ptr<double[]> odata = out()->move_block(c3, c4);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c4);
    sort_indices<0,1,1,1,1,1>(i0data, odata, c3.size(), c4.size());
  }
  out()->put_block(odata, c3, c4);
}

void Task254::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2, c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, c1, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4, c1, x1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), c1.size(), x1.size());
    // tensor label: I47
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a2)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a2.size());
    dgemm_("T", "N", c3.size()*c1.size()*x1.size(), a2.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, c3.size()*c1.size()*x1.size());
  }
  sort_indices<1,3,0,2,1,1,1,1>(odata_sorted, odata, c3.size(), c1.size(), x1.size(), a2.size());
  out()->put_block(odata, c1, a2, c3, x1);
}

void Task255::Task_local::compute() {
  const Index a4 = b(0);
  const Index a2 = b(1);
  // tensor label: I47
  std::unique_ptr<double[]> odata = out()->move_block(a4, a2);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a2);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, a4.size(), a2.size());
  }
  out()->put_block(odata, a4, a2);
}

void Task256::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2, c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c4, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c4, x1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c4.size(), x1.size());
    // tensor label: I50
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, c4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, c4)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), c4.size());
    dgemm_("T", "N", c1.size()*a2.size()*x1.size(), c3.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c1.size()*a2.size()*x1.size());
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), c3.size());
  out()->put_block(odata, c1, a2, c3, x1);
}

void Task257::Task_local::compute() {
  const Index c3 = b(0);
  const Index c4 = b(1);
  // tensor label: I50
  std::unique_ptr<double[]> odata = out()->move_block(c3, c4);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c4);
    sort_indices<0,1,1,1,-2,1>(i0data, odata, c3.size(), c4.size());
  }
  out()->put_block(odata, c3, c4);
}

void Task258::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2, c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, x1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x1.size());
    // tensor label: I53
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a2)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a2.size());
    dgemm_("T", "N", c1.size()*c3.size()*x1.size(), a2.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, c1.size()*c3.size()*x1.size());
  }
  sort_indices<0,3,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), x1.size(), a2.size());
  out()->put_block(odata, c1, a2, c3, x1);
}

void Task259::Task_local::compute() {
  const Index a4 = b(0);
  const Index a2 = b(1);
  // tensor label: I53
  std::unique_ptr<double[]> odata = out()->move_block(a4, a2);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a2);
    sort_indices<0,1,1,1,2,1>(i0data, odata, a4.size(), a2.size());
  }
  out()->put_block(odata, a4, a2);
}

void Task260::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2, c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c5 : *range_[0]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c5, a2, c4, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c5, a2, c4, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c5.size(), a2.size(), c4.size(), x1.size());
      // tensor label: I508
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c5, c3, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c5, c3, c4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), c5.size(), c3.size(), c4.size());
      dgemm_("T", "N", a2.size()*x1.size(), c1.size()*c3.size(), c5.size()*c4.size(),
             1.0, i0data_sorted, c5.size()*c4.size(), i1data_sorted, c5.size()*c4.size(),
             1.0, odata_sorted, a2.size()*x1.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, a2.size(), x1.size(), c1.size(), c3.size());
  out()->put_block(odata, c1, a2, c3, x1);
}

void Task261::Task_local::compute() {
  const Index c1 = b(0);
  const Index c5 = b(1);
  const Index c3 = b(2);
  const Index c4 = b(3);
  // tensor label: I508
  std::unique_ptr<double[]> odata = out()->move_block(c1, c5, c3, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c5, c3, c4);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, c1.size(), c5.size(), c3.size(), c4.size());
  }
  out()->put_block(odata, c1, c5, c3, c4);
}

void Task262::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2, c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& c5 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c5, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a2, c5, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c5.size(), x1.size());
      // tensor label: I511
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c5, c3, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c5, c3, c4)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), c5.size(), c3.size(), c4.size());
      dgemm_("T", "N", a2.size()*x1.size(), c1.size()*c3.size(), c5.size()*c4.size(),
             1.0, i0data_sorted, c5.size()*c4.size(), i1data_sorted, c5.size()*c4.size(),
             1.0, odata_sorted, a2.size()*x1.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, a2.size(), x1.size(), c1.size(), c3.size());
  out()->put_block(odata, c1, a2, c3, x1);
}

void Task263::Task_local::compute() {
  const Index c1 = b(0);
  const Index c5 = b(1);
  const Index c3 = b(2);
  const Index c4 = b(3);
  // tensor label: I511
  std::unique_ptr<double[]> odata = out()->move_block(c1, c5, c3, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c5, c3, c4);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c1.size(), c5.size(), c3.size(), c4.size());
  }
  out()->put_block(odata, c1, c5, c3, c4);
}

void Task264::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2, c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c5 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c5, a4, c3, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c5, a4, c3, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c5.size(), a4.size(), c3.size(), x1.size());
      // tensor label: I514
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c5, a4, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c5, a4, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), c5.size(), a4.size(), a2.size());
      dgemm_("T", "N", c3.size()*x1.size(), c1.size()*a2.size(), c5.size()*a4.size(),
             1.0, i0data_sorted, c5.size()*a4.size(), i1data_sorted, c5.size()*a4.size(),
             1.0, odata_sorted, c3.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c3.size(), x1.size(), c1.size(), a2.size());
  out()->put_block(odata, c1, a2, c3, x1);
}

void Task265::Task_local::compute() {
  const Index c1 = b(0);
  const Index c5 = b(1);
  const Index a4 = b(2);
  const Index a2 = b(3);
  // tensor label: I514
  std::unique_ptr<double[]> odata = out()->move_block(c1, c5, a4, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c5, a4, a2);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c1.size(), c5.size(), a4.size(), a2.size());
  }
  out()->put_block(odata, c1, c5, a4, a2);
}

void Task266::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2, c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a5, c3, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a5, c3, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a5.size(), c3.size(), x1.size());
      // tensor label: I517
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, a5, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, a5, c4)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), a5.size(), c4.size());
      dgemm_("T", "N", c3.size()*x1.size(), c1.size()*a2.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, c3.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c3.size(), x1.size(), c1.size(), a2.size());
  out()->put_block(odata, c1, a2, c3, x1);
}

void Task267::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index a5 = b(2);
  const Index c4 = b(3);
  // tensor label: I517
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2, a5, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, a5, c4);
    sort_indices<0,1,2,3,1,1,4,1>(i0data, odata, c1.size(), a2.size(), a5.size(), c4.size());
  }
  out()->put_block(odata, c1, a2, a5, c4);
}

void Task268::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2, c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c5 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, c5, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4, c5, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), c5.size(), x1.size());
      // tensor label: I520
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c5, a4, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c5, a4, a2)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), c5.size(), a4.size(), a2.size());
      dgemm_("T", "N", c3.size()*x1.size(), c1.size()*a2.size(), c5.size()*a4.size(),
             1.0, i0data_sorted, c5.size()*a4.size(), i1data_sorted, c5.size()*a4.size(),
             1.0, odata_sorted, c3.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c3.size(), x1.size(), c1.size(), a2.size());
  out()->put_block(odata, c1, a2, c3, x1);
}

void Task269::Task_local::compute() {
  const Index c1 = b(0);
  const Index c5 = b(1);
  const Index a4 = b(2);
  const Index a2 = b(3);
  // tensor label: I520
  std::unique_ptr<double[]> odata = out()->move_block(c1, c5, a4, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c5, a4, a2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c1.size(), c5.size(), a4.size(), a2.size());
  }
  out()->put_block(odata, c1, c5, a4, a2);
}

void Task270::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2, c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a5 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a5, c4, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a5, c4, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a5.size(), c4.size(), x1.size());
      // tensor label: I523
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, a5, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, a5, c4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), a5.size(), c4.size());
      dgemm_("T", "N", c3.size()*x1.size(), c1.size()*a2.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, c3.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c3.size(), x1.size(), c1.size(), a2.size());
  out()->put_block(odata, c1, a2, c3, x1);
}

void Task271::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index a5 = b(2);
  const Index c4 = b(3);
  // tensor label: I523
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2, a5, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, a5, c4);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c1.size(), a2.size(), a5.size(), c4.size());
  }
  out()->put_block(odata, c1, a2, a5, c4);
}

void Task272::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2, c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c5 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c5, a4, c1, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c5, a4, c1, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c5.size(), a4.size(), c1.size(), x1.size());
      // tensor label: I526
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, c5, a4, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, c5, a4, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), c5.size(), a4.size(), a2.size());
      dgemm_("T", "N", c1.size()*x1.size(), c3.size()*a2.size(), c5.size()*a4.size(),
             1.0, i0data_sorted, c5.size()*a4.size(), i1data_sorted, c5.size()*a4.size(),
             1.0, odata_sorted, c1.size()*x1.size());
    }
  }
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), x1.size(), c3.size(), a2.size());
  out()->put_block(odata, c1, a2, c3, x1);
}

void Task273::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  const Index a4 = b(2);
  const Index a2 = b(3);
  // tensor label: I526
  std::unique_ptr<double[]> odata = out()->move_block(c3, c5, a4, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c5, a4, a2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c3.size(), c5.size(), a4.size(), a2.size());
  }
  out()->put_block(odata, c3, c5, a4, a2);
}

void Task274::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2, c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a5, c1, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a5, c1, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a5.size(), c1.size(), x1.size());
      // tensor label: I529
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, a5, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, a5, c4)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), a5.size(), c4.size());
      dgemm_("T", "N", c1.size()*x1.size(), c3.size()*a2.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, c1.size()*x1.size());
    }
  }
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), x1.size(), c3.size(), a2.size());
  out()->put_block(odata, c1, a2, c3, x1);
}

void Task275::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index a5 = b(2);
  const Index c4 = b(3);
  // tensor label: I529
  std::unique_ptr<double[]> odata = out()->move_block(c3, a2, a5, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, a5, c4);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c3.size(), a2.size(), a5.size(), c4.size());
  }
  out()->put_block(odata, c3, a2, a5, c4);
}

void Task276::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2, c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c5 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c5, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c5, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c5.size(), x1.size());
      // tensor label: I532
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, c5, a4, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, c5, a4, a2)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), c5.size(), a4.size(), a2.size());
      dgemm_("T", "N", c1.size()*x1.size(), c3.size()*a2.size(), c5.size()*a4.size(),
             1.0, i0data_sorted, c5.size()*a4.size(), i1data_sorted, c5.size()*a4.size(),
             1.0, odata_sorted, c1.size()*x1.size());
    }
  }
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), x1.size(), c3.size(), a2.size());
  out()->put_block(odata, c1, a2, c3, x1);
}

void Task277::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  const Index a4 = b(2);
  const Index a2 = b(3);
  // tensor label: I532
  std::unique_ptr<double[]> odata = out()->move_block(c3, c5, a4, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c5, a4, a2);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c3.size(), c5.size(), a4.size(), a2.size());
  }
  out()->put_block(odata, c3, c5, a4, a2);
}

void Task278::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2, c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a5 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a5, c4, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a5, c4, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a5.size(), c4.size(), x1.size());
      // tensor label: I535
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, a5, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, a5, c4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), a5.size(), c4.size());
      dgemm_("T", "N", c1.size()*x1.size(), c3.size()*a2.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, c1.size()*x1.size());
    }
  }
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), x1.size(), c3.size(), a2.size());
  out()->put_block(odata, c1, a2, c3, x1);
}

void Task279::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index a5 = b(2);
  const Index c4 = b(3);
  // tensor label: I535
  std::unique_ptr<double[]> odata = out()->move_block(c3, a2, a5, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, a5, c4);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c3.size(), a2.size(), a5.size(), c4.size());
  }
  out()->put_block(odata, c3, a2, a5, c4);
}

void Task280::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2, c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a5 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a5, c4, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a5, c4, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a5.size(), c4.size(), a2.size());
      // tensor label: I613
      std::unique_ptr<double[]> i1data = in(1)->get_block(a5, x1, c1, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, x1, c1, c4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, a5.size(), x1.size(), c1.size(), c4.size());
      dgemm_("T", "N", c3.size()*a2.size(), x1.size()*c1.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, c3.size()*a2.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), x1.size(), c1.size());
  out()->put_block(odata, c1, a2, c3, x1);
}

void Task281::Task_local::compute() {
  const Index a5 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index c4 = b(3);
  // tensor label: I613
  std::unique_ptr<double[]> odata = out()->move_block(a5, x1, c1, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, x1, c1, c4);
    sort_indices<0,1,2,3,1,1,-4,1>(i0data, odata, a5.size(), x1.size(), c1.size(), c4.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, x1, a5, c4);
    sort_indices<2,1,0,3,1,1,2,1>(i1data, odata, c1.size(), x1.size(), a5.size(), c4.size());
  }
  out()->put_block(odata, a5, x1, c1, c4);
}

void Task282::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2, c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c4, a5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c4, a5)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c4.size(), a5.size());
      // tensor label: I616
      std::unique_ptr<double[]> i1data = in(1)->get_block(a5, x1, c1, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, x1, c1, c4)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a5.size(), x1.size(), c1.size(), c4.size());
      dgemm_("T", "N", c3.size()*a2.size(), x1.size()*c1.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, c3.size()*a2.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), x1.size(), c1.size());
  out()->put_block(odata, c1, a2, c3, x1);
}

void Task283::Task_local::compute() {
  const Index a5 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index c4 = b(3);
  // tensor label: I616
  std::unique_ptr<double[]> odata = out()->move_block(a5, x1, c1, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, x1, c1, c4);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, a5.size(), x1.size(), c1.size(), c4.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, x1, a5, c4);
    sort_indices<2,1,0,3,1,1,-4,1>(i1data, odata, c1.size(), x1.size(), a5.size(), c4.size());
  }
  out()->put_block(odata, a5, x1, c1, c4);
}

void Task284::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2, c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a5 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a5, c4, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a5, c4, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a5.size(), c4.size(), a2.size());
      // tensor label: I625
      std::unique_ptr<double[]> i1data = in(1)->get_block(a5, x1, c3, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, x1, c3, c4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, a5.size(), x1.size(), c3.size(), c4.size());
      dgemm_("T", "N", c1.size()*a2.size(), x1.size()*c3.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), c3.size());
  out()->put_block(odata, c1, a2, c3, x1);
}

void Task285::Task_local::compute() {
  const Index a5 = b(0);
  const Index x1 = b(1);
  const Index c3 = b(2);
  const Index c4 = b(3);
  // tensor label: I625
  std::unique_ptr<double[]> odata = out()->move_block(a5, x1, c3, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, x1, c3, c4);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, a5.size(), x1.size(), c3.size(), c4.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c3, x1, a5, c4);
    sort_indices<2,1,0,3,1,1,-4,1>(i1data, odata, c3.size(), x1.size(), a5.size(), c4.size());
  }
  out()->put_block(odata, a5, x1, c3, c4);
}

void Task286::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2, c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c4, a5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c4, a5)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c4.size(), a5.size());
      // tensor label: I628
      std::unique_ptr<double[]> i1data = in(1)->get_block(a5, x1, c3, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, x1, c3, c4)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a5.size(), x1.size(), c3.size(), c4.size());
      dgemm_("T", "N", c1.size()*a2.size(), x1.size()*c3.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), c3.size());
  out()->put_block(odata, c1, a2, c3, x1);
}

void Task287::Task_local::compute() {
  const Index a5 = b(0);
  const Index x1 = b(1);
  const Index c3 = b(2);
  const Index c4 = b(3);
  // tensor label: I628
  std::unique_ptr<double[]> odata = out()->move_block(a5, x1, c3, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, x1, c3, c4);
    sort_indices<0,1,2,3,1,1,-4,1>(i0data, odata, a5.size(), x1.size(), c3.size(), c4.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c3, x1, a5, c4);
    sort_indices<2,1,0,3,1,1,8,1>(i1data, odata, c3.size(), x1.size(), a5.size(), c4.size());
  }
  out()->put_block(odata, a5, x1, c3, c4);
}

void Task288::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2, c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a5 : *range_[2]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a5, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a5, c3, a4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a5.size(), c3.size(), a4.size());
      // tensor label: I637
      std::unique_ptr<double[]> i1data = in(1)->get_block(a5, x1, a4, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, x1, a4, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, a5.size(), x1.size(), a4.size(), a2.size());
      dgemm_("T", "N", c1.size()*c3.size(), x1.size()*a2.size(), a5.size()*a4.size(),
             1.0, i0data_sorted, a5.size()*a4.size(), i1data_sorted, a5.size()*a4.size(),
             1.0, odata_sorted, c1.size()*c3.size());
    }
  }
  sort_indices<0,3,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), x1.size(), a2.size());
  out()->put_block(odata, c1, a2, c3, x1);
}

void Task289::Task_local::compute() {
  const Index a5 = b(0);
  const Index x1 = b(1);
  const Index a4 = b(2);
  const Index a2 = b(3);
  // tensor label: I637
  std::unique_ptr<double[]> odata = out()->move_block(a5, x1, a4, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, x1, a4, a2);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, a5.size(), x1.size(), a4.size(), a2.size());
  }
  out()->put_block(odata, a5, x1, a4, a2);
}

void Task290::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2, c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a5)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a5.size());
      // tensor label: I640
      std::unique_ptr<double[]> i1data = in(1)->get_block(a5, x1, a4, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, x1, a4, a2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a5.size(), x1.size(), a4.size(), a2.size());
      dgemm_("T", "N", c1.size()*c3.size(), x1.size()*a2.size(), a5.size()*a4.size(),
             1.0, i0data_sorted, a5.size()*a4.size(), i1data_sorted, a5.size()*a4.size(),
             1.0, odata_sorted, c1.size()*c3.size());
    }
  }
  sort_indices<0,3,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), x1.size(), a2.size());
  out()->put_block(odata, c1, a2, c3, x1);
}

void Task291::Task_local::compute() {
  const Index a5 = b(0);
  const Index x1 = b(1);
  const Index a4 = b(2);
  const Index a2 = b(3);
  // tensor label: I640
  std::unique_ptr<double[]> odata = out()->move_block(a5, x1, a4, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, x1, a4, a2);
    sort_indices<0,1,2,3,1,1,4,1>(i0data, odata, a5.size(), x1.size(), a4.size(), a2.size());
  }
  out()->put_block(odata, a5, x1, a4, a2);
}

void Task292::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), x1.size());
    // tensor label: I55
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c3, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c3, x1, x0)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), c3.size(), x1.size(), x0.size());
    dgemm_("T", "N", c1.size(), a2.size()*c3.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size());
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task293::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I55
  std::unique_ptr<double[]> odata = out()->move_block(a2, c3, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma18
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x0, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x0.size(), x2.size());
      // tensor label: I56
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

void Task294::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x2 = b(3);
  // tensor label: I56
  std::unique_ptr<double[]> odata = out()->move_block(x3, a2, c3, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c3, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), a2.size(), c3.size(), x2.size());
  }
  out()->put_block(odata, x3, a2, c3, x2);
}

void Task295::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I55
  std::unique_ptr<double[]> odata = out()->move_block(a2, c3, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma10
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: I62
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

void Task296::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I62
  std::unique_ptr<double[]> odata = out()->move_block(c3, a2, x3, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, x3, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c3.size(), a2.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, c3, a2, x3, x2);
}

void Task297::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), x1.size());
    // tensor label: I58
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c1, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c1, x0, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a2.size(), c1.size(), x0.size(), x1.size());
    dgemm_("T", "N", c3.size(), a2.size()*c1.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c3.size());
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task298::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I58
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma10
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: I59
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

void Task299::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x2 = b(3);
  // tensor label: I59
  std::unique_ptr<double[]> odata = out()->move_block(x3, a2, c1, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), a2.size(), c1.size(), x2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, x3, x2);
    sort_indices<2,1,0,3,1,1,2,1>(i1data, odata, c1.size(), a2.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, x3, a2, c1, x2);
}

#endif
