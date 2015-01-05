//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks21.cc
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


#include <src/smith/CASPT2_tasks21.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task1000::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1164
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
    dscal_(x1.size()*x0.size()*a2.size()*c1.size(), e0_, i0data.get(), 1);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, c1.size(), a2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x1, x0, a2, c1);
}

void Task1001::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I917
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x3, x2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x3.size(), x2.size());
      // tensor label: I1167
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, a2, c1)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), a2.size(), c1.size());
      dgemm_("T", "N", x3.size()*x2.size(), x1.size()*x0.size(), a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, x3, x2);
}

void Task1002::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1167
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
    dscal_(x1.size()*x0.size()*a2.size()*c1.size(), e0_, i0data.get(), 1);
    sort_indices<3,2,1,0,1,1,-1,1>(i0data, odata, c1.size(), a2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x1, x0, a2, c1);
}

void Task1003::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I917
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x3, x2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x3.size(), x2.size());
      // tensor label: I1197
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c2, a1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c2, a1, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), c2.size(), a1.size(), x0.size());
      dgemm_("T", "N", x3.size()*x2.size(), x1.size()*x0.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, x3, x2);
}

void Task1004::Task_local::compute() {
  const Index x1 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1197
  std::unique_ptr<double[]> odata = out()->move_block(x1, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, x0.size(), a1.size(), c2.size(), x1.size());
  }
  out()->put_block(odata, x1, c2, a1, x0);
}

void Task1005::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I917
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, c2, a1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), c2.size(), a1.size());
      // tensor label: I1206
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c2, a1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c2, a1, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), c2.size(), a1.size(), x0.size());
      dgemm_("T", "N", x3.size()*x2.size(), x1.size()*x0.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, x3, x2);
}

void Task1006::Task_local::compute() {
  const Index x1 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1206
  std::unique_ptr<double[]> odata = out()->move_block(x1, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, x0.size(), a1.size(), c2.size(), x1.size());
  }
  out()->put_block(odata, x1, c2, a1, x0);
}

void Task1007::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I917
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x3, x2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x3.size(), x2.size());
      // tensor label: I1209
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, a2, c1)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), a2.size(), c1.size());
      dgemm_("T", "N", x3.size()*x2.size(), x1.size()*x0.size(), a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, x3, x2);
}

void Task1008::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1209
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, c1.size(), a2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x1, x0, a2, c1);
}

void Task1009::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I917
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c1, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c1.size(), x2.size());
      // tensor label: I1215
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, a2, c1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), a2.size(), c1.size());
      dgemm_("T", "N", x3.size()*x2.size(), x1.size()*x0.size(), a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, x3, x2);
}

void Task1010::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1215
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, c1.size(), a2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x1, x0, a2, c1);
}

void Task1011::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I917
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, c1, a2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), c1.size(), a2.size());
      // tensor label: I1218
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, a2, c1)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), a2.size(), c1.size());
      dgemm_("T", "N", x3.size()*x2.size(), x1.size()*x0.size(), a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, x3, x2);
}

void Task1012::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1218
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, c1.size(), a2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x1, x0, a2, c1);
}

void Task1013::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I917
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x0.size(), x1.size());
      // tensor label: I1251
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c2, a1, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c2, a1, x3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), c2.size(), a1.size(), x3.size());
      dgemm_("T", "N", x0.size()*x1.size(), x2.size()*x3.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size());
  out()->put_block(odata, x1, x0, x3, x2);
}

void Task1014::Task_local::compute() {
  const Index x2 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x3 = b(3);
  // tensor label: I1251
  std::unique_ptr<double[]> odata = out()->move_block(x2, c2, a1, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, x2);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, x3.size(), a1.size(), c2.size(), x2.size());
  }
  out()->put_block(odata, x2, c2, a1, x3);
}

void Task1015::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I917
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1, c2, a1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size(), c2.size(), a1.size());
      // tensor label: I1260
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c2, a1, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c2, a1, x3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), c2.size(), a1.size(), x3.size());
      dgemm_("T", "N", x0.size()*x1.size(), x2.size()*x3.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size());
  out()->put_block(odata, x1, x0, x3, x2);
}

void Task1016::Task_local::compute() {
  const Index x2 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x3 = b(3);
  // tensor label: I1260
  std::unique_ptr<double[]> odata = out()->move_block(x2, c2, a1, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, x2);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, x3.size(), a1.size(), c2.size(), x2.size());
  }
  out()->put_block(odata, x2, c2, a1, x3);
}

void Task1017::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I917
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x0.size(), x1.size());
      // tensor label: I1263
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x3, a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x3, a2, c1)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, x2.size(), x3.size(), a2.size(), c1.size());
      dgemm_("T", "N", x0.size()*x1.size(), x2.size()*x3.size(), a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size());
  out()->put_block(odata, x1, x0, x3, x2);
}

void Task1018::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1263
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x3, x2);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, c1.size(), a2.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, x2, x3, a2, c1);
}

void Task1019::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I917
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a2, c1, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a2, c1, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a2.size(), c1.size(), x1.size());
      // tensor label: I1269
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x3, a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x3, a2, c1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, x2.size(), x3.size(), a2.size(), c1.size());
      dgemm_("T", "N", x0.size()*x1.size(), x2.size()*x3.size(), a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size());
  out()->put_block(odata, x1, x0, x3, x2);
}

void Task1020::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1269
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x3, x2);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, c1.size(), a2.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, x2, x3, a2, c1);
}

void Task1021::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I917
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1, c1, a2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size(), c1.size(), a2.size());
      // tensor label: I1272
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x3, a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x3, a2, c1)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, x2.size(), x3.size(), a2.size(), c1.size());
      dgemm_("T", "N", x0.size()*x1.size(), x2.size()*x3.size(), a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size());
  out()->put_block(odata, x1, x0, x3, x2);
}

void Task1022::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1272
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x3, x2);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, c1.size(), a2.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, x2, x3, a2, c1);
}

void Task1023::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I782
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        for (auto& x3 : *range_[1]) {
          for (auto& x1 : *range_[1]) {
            for (auto& x2 : *range_[1]) {
              // tensor label: Gamma309
              std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x3, x1, x2);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x5, x0, x4, x3, x1, x2)]);
              sort_indices<1,2,3,4,5,6,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x1.size(), x2.size());
              // tensor label: I925
              std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, x2, x5, x4, x3);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, x2, x5, x4, x3)]);
              sort_indices<3,1,4,5,0,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), x2.size(), x5.size(), x4.size(), x3.size());
              dgemm_("T", "N", ci0.size(), 1, x1.size()*x0.size()*x2.size()*x5.size()*x4.size()*x3.size(),
                     1.0, i0data_sorted, x1.size()*x0.size()*x2.size()*x5.size()*x4.size()*x3.size(), i1data_sorted, x1.size()*x0.size()*x2.size()*x5.size()*x4.size()*x3.size(),
                     1.0, odata_sorted, ci0.size());
            }
          }
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1024::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x5 = b(3);
  const Index x4 = b(4);
  const Index x3 = b(5);
  // tensor label: I925
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x2, x5, x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x2, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x2, x5, x4, x3), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, x3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
    // tensor label: I926
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a1, x0, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a1, x0, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), a1.size(), x0.size(), x2.size());
    dgemm_("T", "N", x5.size()*x4.size()*x3.size(), x1.size()*x0.size()*x2.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x5.size()*x4.size()*x3.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), x1.size(), x0.size(), x2.size());
  out()->put_block(odata, x1, x0, x2, x5, x4, x3);
}

void Task1025::Task_local::compute() {
  const Index x1 = b(0);
  const Index a1 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  // tensor label: I926
  std::unique_ptr<double[]> odata = out()->move_block(x1, a1, x0, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, a1, x0, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, a1, x0, x2), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), x2.size());
    // tensor label: I927
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c2, a1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c2, a1, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), c2.size(), a1.size(), x0.size());
    dgemm_("T", "N", x2.size(), x1.size()*a1.size()*x0.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, x2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), a1.size(), x0.size());
  out()->put_block(odata, x1, a1, x0, x2);
}

void Task1026::Task_local::compute() {
  const Index x1 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I927
  std::unique_ptr<double[]> odata = out()->move_block(x1, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, x0.size(), a1.size(), c2.size(), x1.size());
  }
  out()->put_block(odata, x1, c2, a1, x0);
}

void Task1027::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I782
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma310
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x1, x0)]);
      sort_indices<1,2,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x1.size(), x0.size());
      // tensor label: I929
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
      dgemm_("T", "N", ci0.size(), 1, x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, ci0.size());
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1028::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I929
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, x1)]);
      sort_indices<2,1,0,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), x1.size());
      // tensor label: I930
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size());
      dgemm_("T", "N", x1.size()*x0.size(), 1, c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->put_block(odata, x1, x0);
}

void Task1029::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  // tensor label: I930
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a4, c3, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a4, c3, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a4.size(), c3.size(), a1.size());
      // tensor label: I931
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

void Task1030::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  // tensor label: I931
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3);
    sort_indices<0,1,1,1,1,1>(i0data, odata, a4.size(), c3.size());
  }
  out()->put_block(odata, a4, c3);
}

void Task1031::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  // tensor label: I930
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c3, a4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c3.size(), a4.size());
      // tensor label: I935
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

void Task1032::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  // tensor label: I935
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3);
    sort_indices<0,1,1,1,-2,1>(i0data, odata, a4.size(), c3.size());
  }
  out()->put_block(odata, a4, c3);
}

void Task1033::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I929
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x0, x1)]);
      sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x0.size(), x1.size());
      // tensor label: I984
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size());
      dgemm_("T", "N", x1.size()*x0.size(), 1, c1.size()*a2.size(),
             1.0, i0data_sorted, c1.size()*a2.size(), i1data_sorted, c1.size()*a2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->put_block(odata, x1, x0);
}

void Task1034::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: I984
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
      // tensor label: I985
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size());
      dgemm_("T", "N", c1.size()*a2.size(), 1, a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size());
  out()->put_block(odata, c1, a2);
}

void Task1035::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  // tensor label: I985
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3);
    sort_indices<0,1,1,1,-2,1>(i0data, odata, a4.size(), c3.size());
  }
  out()->put_block(odata, a4, c3);
}

void Task1036::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: I984
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
      // tensor label: I989
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size());
      dgemm_("T", "N", c1.size()*a2.size(), 1, a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size());
  out()->put_block(odata, c1, a2);
}

void Task1037::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  // tensor label: I989
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3);
    sort_indices<0,1,1,1,4,1>(i0data, odata, a4.size(), c3.size());
  }
  out()->put_block(odata, a4, c3);
}

void Task1038::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I929
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c1, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c1.size(), x0.size());
      // tensor label: I1035
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c1)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), c1.size());
      dgemm_("T", "N", x1.size()*x0.size(), 1, a4.size()*c1.size(),
             1.0, i0data_sorted, a4.size()*c1.size(), i1data_sorted, a4.size()*c1.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x1, x0);
}

void Task1039::Task_local::compute() {
  const Index a4 = b(0);
  const Index c1 = b(1);
  // tensor label: I1035
  std::unique_ptr<double[]> odata = out()->move_block(a4, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size());
      // tensor label: I1036
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3, a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3, a2, c1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size(), a2.size(), c1.size());
      dgemm_("T", "N", 1, a4.size()*c1.size(), c3.size()*a2.size(),
             1.0, i0data_sorted, c3.size()*a2.size(), i1data_sorted, c3.size()*a2.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size());
  out()->put_block(odata, a4, c1);
}

void Task1040::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1036
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->put_block(odata, a4, c3, a2, c1);
}

void Task1041::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I929
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2, c1, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), x0.size());
      // tensor label: I1039
      std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c1)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a2.size(), c1.size());
      dgemm_("T", "N", x1.size()*x0.size(), 1, a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x1, x0);
}

void Task1042::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  // tensor label: I1039
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size());
      // tensor label: I1040
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3, a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3, a2, c1)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size(), a2.size(), c1.size());
      dgemm_("T", "N", 1, a2.size()*c1.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size());
  out()->put_block(odata, a2, c1);
}

void Task1043::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1040
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<3,2,1,0,1,1,-2,1>(i0data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->put_block(odata, a4, c3, a2, c1);
}

void Task1044::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I929
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), x1.size(), x0.size());
      // tensor label: I1043
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c1)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a4.size(), c1.size());
      dgemm_("T", "N", x1.size()*x0.size(), 1, a4.size()*c1.size(),
             1.0, i0data_sorted, a4.size()*c1.size(), i1data_sorted, a4.size()*c1.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x1, x0);
}

void Task1045::Task_local::compute() {
  const Index a4 = b(0);
  const Index c1 = b(1);
  // tensor label: I1043
  std::unique_ptr<double[]> odata = out()->move_block(a4, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size());
      // tensor label: I1044
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3, a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3, a2, c1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size(), a2.size(), c1.size());
      dgemm_("T", "N", 1, a4.size()*c1.size(), c3.size()*a2.size(),
             1.0, i0data_sorted, c3.size()*a2.size(), i1data_sorted, c3.size()*a2.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size());
  out()->put_block(odata, a4, c1);
}

void Task1046::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1044
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<3,2,1,0,1,1,-2,1>(i0data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->put_block(odata, a4, c3, a2, c1);
}

void Task1047::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I929
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x1.size(), x0.size());
      // tensor label: I1047
      std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c1)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a2.size(), c1.size());
      dgemm_("T", "N", x1.size()*x0.size(), 1, a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x1, x0);
}

void Task1048::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  // tensor label: I1047
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size());
      // tensor label: I1048
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3, a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3, a2, c1)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size(), a2.size(), c1.size());
      dgemm_("T", "N", 1, a2.size()*c1.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size());
  out()->put_block(odata, a2, c1);
}

void Task1049::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1048
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<3,2,1,0,1,1,4,1>(i0data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->put_block(odata, a4, c3, a2, c1);
}

