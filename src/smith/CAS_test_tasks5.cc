//
// BAGEL - Parallel electron correlation program.
// Filename: CAS_test_tasks5.cc
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


#include <src/smith/CAS_test_tasks5.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CAS_test;

void Task200::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I116
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task201::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a3 = b(3);
  const Index c2 = b(4);
  const Index a1 = b(5);
  // tensor label: I115
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, a3, c2, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());
    // tensor label: I119
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x1, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x1, x2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), x0.size()*x1.size()*x2.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }
  sort_indices<5,4,3,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size(), x1.size(), x2.size());
  out()->put_block(odata, x2, x1, x0, a3, c2, a1);
}

void Task202::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: I119
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x1, x2);
  {
    // tensor label: Gamma32
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x1, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x0.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, x3, x0, x1, x2);
}

void Task203::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I80
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), x2.size());
      // tensor label: I285
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      dgemm_("T", "N", a1.size()*c2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), x0.size(), x1.size());
  out()->put_block(odata, x1, x0, c2, a1);
}

void Task204::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: I285
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x1, x2);
  {
    // tensor label: Gamma32
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x1, x2);
    dscal_(x3.size()*x0.size()*x1.size()*x2.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), x0.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, x3, x0, x1, x2);
}

void Task205::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I80
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x3.size(), x2.size());
      // tensor label: I287
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c2, a1);
}

void Task206::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I287
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    dscal_(x3.size()*x2.size()*x1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task207::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I80
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x3.size(), x2.size());
      // tensor label: I313
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c2, a1);
}

void Task208::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I313
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task209::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I80
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x3, x2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x3, x2, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), x3.size(), x2.size(), a1.size());
      // tensor label: I315
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x3, x2, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x3, x2, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), x3.size(), x2.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c2, a1);
}

void Task210::Task_local::compute() {
  const Index x1 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x0 = b(3);
  // tensor label: I315
  std::unique_ptr<double[]> odata = out()->move_block(x1, x3, x2, x0);
  {
    // tensor label: Gamma29
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x3, x2, x0);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x1.size(), x3.size(), x2.size(), x0.size());
  }
  out()->put_block(odata, x1, x3, x2, x0);
}

void Task211::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I80
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), x2.size());
      // tensor label: I317
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      dgemm_("T", "N", a1.size()*c2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), x0.size(), x1.size());
  out()->put_block(odata, x1, x0, c2, a1);
}

void Task212::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: I317
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x1, x2);
  {
    // tensor label: Gamma32
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x1, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x0.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, x3, x0, x1, x2);
}

void Task213::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I80
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, c2, a1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), c2.size(), a1.size());
      // tensor label: I319
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c2, a1);
}

void Task214::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I319
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task215::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I80
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a1), 0.0);
  // tensor label: h1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size());
  // tensor label: I343
  std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
  dgemm_("T", "N", c2.size()*a1.size(), x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a1.size());
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c2, a1);
}

void Task216::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I343
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-2,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task217::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, c1, a2);
  {
    // tensor label: I120
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0, c1, a2);
    sort_indices<1,0,2,3,1,1,1,1>(i0data, odata, x1.size(), x0.size(), c1.size(), a2.size());
  }
  out()->put_block(odata, x0, x1, c1, a2);
}

void Task218::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I120
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), a2.size());
    // tensor label: I121
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, x0, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, x0, c1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), x0.size(), c1.size());
    dgemm_("T", "N", a2.size(), x1.size()*x0.size()*c1.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, x1, x0, c1, a2);
}

void Task219::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I121
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, c1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c1, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, c1, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), c1.size(), x3.size());
        // tensor label: I122
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x2, x3, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x2, x3, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
        dgemm_("T", "N", c1.size(), x2.size()*x1.size()*x0.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c1.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x1, x0, c1);
}

void Task220::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x2 = b(2);
  const Index x3 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I122
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, x2, x3, x1, x0);
  {
    // tensor label: Gamma6
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x2, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x5, x4, x2, x3, x1, x0);
}

void Task221::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I120
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, c3)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), c3.size());
      // tensor label: I124
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, x0, c3, a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, x0, c3, a2, c1)]);
      sort_indices<0,3,1,2,4,5,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), x0.size(), c3.size(), a2.size(), c1.size());
      dgemm_("T", "N", 1, x1.size()*x0.size()*a2.size()*c1.size(), x2.size()*c3.size(),
             1.0, i0data_sorted, x2.size()*c3.size(), i1data_sorted, x2.size()*c3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a2.size(), c1.size());
  out()->put_block(odata, x1, x0, c1, a2);
}

void Task222::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);
  const Index c1 = b(5);
  // tensor label: I124
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, c3, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, c3, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, c3, a2, c1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c1, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c1.size(), x3.size());
    // tensor label: I125
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x3, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x3, x1, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
    dgemm_("T", "N", c3.size()*a2.size()*c1.size(), x2.size()*x1.size()*x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c3.size()*a2.size()*c1.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x1, x0, c3, a2, c1);
}

void Task223::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I125
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, x1, x0);
  {
    // tensor label: Gamma7
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x2, x3, x1, x0);
}

void Task224::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);
  const Index c1 = b(5);
  // tensor label: I124
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, c3, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, c3, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, c3, a2, c1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x3.size());
    // tensor label: I128
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x3, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x3, x1, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
    dgemm_("T", "N", c1.size()*a2.size()*c3.size(), x2.size()*x1.size()*x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<3,4,5,2,1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x1, x0, c3, a2, c1);
}

void Task225::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I128
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, x1, x0);
  {
    // tensor label: Gamma7
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x2, x3, x1, x0);
}

void Task226::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I120
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, c1, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, c1, x4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), c1.size(), x4.size());
      // tensor label: I130
      std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x1.size(), x0.size());
      dgemm_("T", "N", a2.size()*c1.size(), x1.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c1, a2);
}

void Task227::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I130
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, x1, x0);
  {
    // tensor label: Gamma34
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x5.size(), x4.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x5, x4, x1, x0);
}

void Task228::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I120
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, c3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), c3.size());
    // tensor label: I132
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, a2, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, a2, c3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), a2.size(), c3.size());
    dgemm_("T", "N", c1.size(), x1.size()*x0.size()*a2.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, c1.size());
  }
  sort_indices<1,2,0,3,1,1,1,1>(odata_sorted, odata, c1.size(), x1.size(), x0.size(), a2.size());
  out()->put_block(odata, x1, x0, c1, a2);
}

void Task229::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c3 = b(3);
  // tensor label: I132
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, a2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a2, c3), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c3, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c3.size(), x2.size());
      // tensor label: I133
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", a2.size()*c3.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a2.size()*c3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a2.size(), c3.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, a2, c3);
}

void Task230::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I133
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task231::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c3 = b(3);
  // tensor label: I132
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, a2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a2, c3), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), x3.size(), x2.size());
      // tensor label: I141
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c3.size()*a2.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c3.size()*a2.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, a2, c3);
}

void Task232::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I141
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task233::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I120
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), a2.size());
    // tensor label: I135
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, a3, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, a3, c1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), a3.size(), c1.size());
    dgemm_("T", "N", a2.size(), x1.size()*x0.size()*c1.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, x1, x0, c1, a2);
}

void Task234::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a3 = b(2);
  const Index c1 = b(3);
  // tensor label: I135
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a3, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, a3, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a3, c1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c1.size(), x2.size());
      // tensor label: I136
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", a3.size()*c1.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a3.size()*c1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a3.size(), c1.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, a3, c1);
}

void Task235::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I136
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task236::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a3 = b(2);
  const Index c1 = b(3);
  // tensor label: I135
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a3, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, a3, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a3, c1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), x3.size(), x2.size());
      // tensor label: I144
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a3.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c1.size()*a3.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a3.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, a3, c1);
}

void Task237::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task238::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I120
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x5, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x5, x4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x5.size(), x4.size());
      // tensor label: I138
      std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), x1.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c1, a2);
}

void Task239::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I138
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, x1, x0);
  {
    // tensor label: Gamma34
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x1, x0);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x5.size(), x4.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x5, x4, x1, x0);
}

void Task240::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I120
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), x2.size());
    // tensor label: I146
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, x0, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, x0, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), x0.size(), a2.size());
    dgemm_("T", "N", c1.size(), x1.size()*x0.size()*a2.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, c1.size());
  }
  sort_indices<1,2,0,3,1,1,1,1>(odata_sorted, odata, c1.size(), x1.size(), x0.size(), a2.size());
  out()->put_block(odata, x1, x0, c1, a2);
}

void Task241::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I146
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, x4, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), x4.size(), x3.size());
        // tensor label: I147
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x2, x4, x3, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x2, x4, x3, x1, x0)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
        dgemm_("T", "N", a2.size(), x2.size()*x1.size()*x0.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x1, x0, a2);
}

void Task242::Task_local::compute() {
  const Index x5 = b(0);
  const Index x2 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I147
  std::unique_ptr<double[]> odata = out()->move_block(x5, x2, x4, x3, x1, x0);
  {
    // tensor label: Gamma51
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2, x4, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x5, x2, x4, x3, x1, x0);
}

void Task243::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I120
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, c3)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), c3.size());
      // tensor label: I149
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, c1, a4, c3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, c1, a4, c3, a2)]);
      sort_indices<3,4,0,1,2,5,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), c1.size(), a4.size(), c3.size(), a2.size());
      dgemm_("T", "N", 1, x1.size()*x0.size()*c1.size()*a2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), a2.size());
  out()->put_block(odata, x1, x0, c1, a2);
}

void Task244::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a4 = b(3);
  const Index c3 = b(4);
  const Index a2 = b(5);
  // tensor label: I149
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a4, c3, a2), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
  // tensor label: I150
  std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
  dgemm_("T", "N", c1.size()*a4.size()*c3.size()*a2.size(), x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a4.size()*c3.size()*a2.size());
  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c1, a4, c3, a2);
}

void Task245::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I150
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-4,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task246::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a4 = b(3);
  const Index c3 = b(4);
  const Index a2 = b(5);
  // tensor label: I149
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a4, c3, a2), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
  // tensor label: I153
  std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
  dgemm_("T", "N", c1.size()*a2.size()*c3.size()*a4.size(), x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a2.size()*c3.size()*a4.size());
  sort_indices<4,5,0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c1, a4, c3, a2);
}

void Task247::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I153
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,8,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task248::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I120
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, x2)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), x2.size());
      // tensor label: I155
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, x0, a3, c1, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, x0, a3, c1, a2)]);
      sort_indices<3,0,1,2,4,5,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), x0.size(), a3.size(), c1.size(), a2.size());
      dgemm_("T", "N", 1, x1.size()*x0.size()*c1.size()*a2.size(), x2.size()*a3.size(),
             1.0, i0data_sorted, x2.size()*a3.size(), i1data_sorted, x2.size()*a3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), a2.size());
  out()->put_block(odata, x1, x0, c1, a2);
}

void Task249::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a3 = b(3);
  const Index c1 = b(4);
  const Index a2 = b(5);
  // tensor label: I155
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, a3, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, a3, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, a3, c1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c1, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c1.size(), a2.size());
    // tensor label: I156
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
    dgemm_("T", "N", a3.size()*c1.size()*a2.size(), x2.size()*x1.size()*x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a3.size()*c1.size()*a2.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c1.size(), a2.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x1, x0, a3, c1, a2);
}

