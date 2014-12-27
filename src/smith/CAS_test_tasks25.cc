//
// BAGEL - Parallel electron correlation program.
// Filename: CAS_test_tasks25.cc
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


#include <src/smith/CAS_test_tasks25.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CAS_test;

void Task1200::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1305
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c2, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c2, x4)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c2.size(), x4.size());
      // tensor label: I1318
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x0, x1, x4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x0, x1, x4)]);
      sort_indices<4,1,0,2,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x0.size(), x1.size(), x4.size());
      dgemm_("T", "N", a1.size()*c2.size(), ci0.size()*x0.size()*x1.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<2,4,3,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), ci0.size(), x0.size(), x1.size());
  out()->put_block(odata, ci0, x1, x0, c2, a1);
}

void Task1201::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  const Index x4 = b(4);
  // tensor label: I1318
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x1, x4);
  {
    // tensor label: Gamma409
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x1, x4);
    sort_indices<0,1,2,3,4,1,1,1,4>(i0data, odata, ci0.size(), x5.size(), x0.size(), x1.size(), x4.size());
  }
  out()->put_block(odata, ci0, x5, x0, x1, x4);
}

void Task1202::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1305
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a1), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, c3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), c3.size());
    // tensor label: I1321
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1, a1, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1, a1, c3)]);
    sort_indices<4,0,1,2,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size(), a1.size(), c3.size());
    dgemm_("T", "N", c2.size(), ci0.size()*x0.size()*x1.size()*a1.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, c2.size());
  }
  sort_indices<1,3,2,0,4,1,1,1,1>(odata_sorted, odata, c2.size(), ci0.size(), x0.size(), x1.size(), a1.size());
  out()->put_block(odata, ci0, x1, x0, c2, a1);
}

void Task1203::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  const Index c3 = b(4);
  // tensor label: I1321
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, c3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c3, x2)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c3.size(), x2.size());
      // tensor label: I1322
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x1, x2)]);
      sort_indices<4,1,0,2,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x1.size(), x2.size());
      dgemm_("T", "N", a1.size()*c3.size(), ci0.size()*x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*c3.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), c3.size(), ci0.size(), x0.size(), x1.size());
  out()->put_block(odata, ci0, x0, x1, a1, c3);
}

void Task1204::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  const Index x2 = b(4);
  // tensor label: I1322
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x1, x2);
  {
    // tensor label: Gamma410
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x1, x2);
    sort_indices<0,1,2,3,4,1,1,-1,4>(i0data, odata, ci0.size(), x3.size(), x0.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, ci0, x3, x0, x1, x2);
}

void Task1205::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  const Index c3 = b(4);
  // tensor label: I1321
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, c3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a1, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a1, x3, x2)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a1.size(), x3.size(), x2.size());
      // tensor label: I1333
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
      sort_indices<2,1,0,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c3.size()*a1.size(), ci0.size()*x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c3.size()*a1.size());
    }
  }
  sort_indices<2,4,3,1,0,1,1,1,1>(odata_sorted, odata, c3.size(), a1.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x0, x1, a1, c3);
}

void Task1206::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1333
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,1,1,1,4>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1207::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1305
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a1), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), a1.size());
    // tensor label: I1325
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1, a3, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1, a3, c2)]);
    sort_indices<3,0,1,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size(), a3.size(), c2.size());
    dgemm_("T", "N", a1.size(), ci0.size()*x0.size()*x1.size()*c2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<1,3,2,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x0.size(), x1.size(), c2.size());
  out()->put_block(odata, ci0, x1, x0, c2, a1);
}

void Task1208::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a3 = b(3);
  const Index c2 = b(4);
  // tensor label: I1325
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a3, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a3, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a3, c2), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, x2)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), x2.size());
      // tensor label: I1326
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x1, x2)]);
      sort_indices<4,1,0,2,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x1.size(), x2.size());
      dgemm_("T", "N", a3.size()*c2.size(), ci0.size()*x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a3.size()*c2.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), ci0.size(), x0.size(), x1.size());
  out()->put_block(odata, ci0, x0, x1, a3, c2);
}

void Task1209::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  const Index x2 = b(4);
  // tensor label: I1326
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x1, x2);
  {
    // tensor label: Gamma410
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x1, x2);
    sort_indices<0,1,2,3,4,1,1,1,4>(i0data, odata, ci0.size(), x3.size(), x0.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, ci0, x3, x0, x1, x2);
}

void Task1210::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a3 = b(3);
  const Index c2 = b(4);
  // tensor label: I1325
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a3, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a3, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a3, c2), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, x3, x2)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), x3.size(), x2.size());
      // tensor label: I1337
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
      sort_indices<2,1,0,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a3.size(), ci0.size()*x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<2,4,3,1,0,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x0, x1, a3, c2);
}

void Task1211::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1337
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,1,1,-1,4>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1212::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1305
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x5, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x5, x4)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x5.size(), x4.size());
      // tensor label: I1329
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x4, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x4, x1, x0)]);
      sort_indices<2,1,0,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x4.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), ci0.size()*x1.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c2, a1);
}

void Task1213::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1329
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x1, x0);
  {
    // tensor label: Gamma412
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x4, x1, x0);
    sort_indices<0,1,2,3,4,1,1,-1,4>(i0data, odata, ci0.size(), x5.size(), x4.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x5, x4, x1, x0);
}

void Task1214::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1305
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), x2.size());
    // tensor label: I1340
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1, x2, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1, x2, a1)]);
    sort_indices<3,0,1,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size(), x2.size(), a1.size());
    dgemm_("T", "N", c2.size(), ci0.size()*x0.size()*x1.size()*a1.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, c2.size());
  }
  sort_indices<1,3,2,0,4,1,1,1,1>(odata_sorted, odata, c2.size(), ci0.size(), x0.size(), x1.size(), a1.size());
  out()->put_block(odata, ci0, x1, x0, c2, a1);
}

void Task1215::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1340
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, x2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, x2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, x2, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, x3)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
        // tensor label: I1341
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x0, x4, x3, x1, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x0, x4, x3, x1, x2)]);
        sort_indices<4,3,1,0,2,5,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x1.size(), x2.size());
        dgemm_("T", "N", a1.size(), ci0.size()*x0.size()*x1.size()*x2.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x0.size(), x1.size(), x2.size());
  out()->put_block(odata, ci0, x0, x1, x2, a1);
}

void Task1216::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x1 = b(5);
  const Index x2 = b(6);
  // tensor label: I1341
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x3, x1, x2);
  {
    // tensor label: Gamma415
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x3, x1, x2);
    sort_indices<0,1,2,3,4,5,6,1,1,1,4>(i0data, odata, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, ci0, x5, x0, x4, x3, x1, x2);
}

void Task1217::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1305
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, c3)]);
      sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a4.size(), c3.size());
      // tensor label: I1344
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0, c2, a4, c3, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0, c2, a4, c3, a1)]);
      sort_indices<5,4,0,1,2,3,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size(), c2.size(), a4.size(), c3.size(), a1.size());
      dgemm_("T", "N", 1, ci0.size()*x1.size()*x0.size()*c2.size()*a1.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,3,4,1,1,1,1>(odata_sorted, odata, ci0.size(), x1.size(), x0.size(), c2.size(), a1.size());
  out()->put_block(odata, ci0, x1, x0, c2, a1);
}

void Task1218::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a4 = b(4);
  const Index c3 = b(5);
  const Index a1 = b(6);
  // tensor label: I1344
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a4, c3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a4, c3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a4, c3, a1), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a4, c3, a1);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a4, c3, a1)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a4.size(), c3.size(), a1.size());
  // tensor label: I1345
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
  sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
  dgemm_("T", "N", c2.size()*a4.size()*c3.size()*a1.size(), ci0.size()*x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a4.size()*c3.size()*a1.size());
  sort_indices<4,5,6,0,1,2,3,1,1,1,1>(odata_sorted, odata, c2.size(), a4.size(), c3.size(), a1.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c2, a4, c3, a1);
}

void Task1219::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1345
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1220::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a4 = b(4);
  const Index c3 = b(5);
  const Index a1 = b(6);
  // tensor label: I1344
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a4, c3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a4, c3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a4, c3, a1), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c3, a4)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c3.size(), a4.size());
  // tensor label: I1349
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
  sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
  dgemm_("T", "N", c2.size()*a1.size()*c3.size()*a4.size(), ci0.size()*x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a1.size()*c3.size()*a4.size());
  sort_indices<4,5,6,0,3,2,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), c3.size(), a4.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c2, a4, c3, a1);
}

void Task1221::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1349
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1222::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1305
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, x2)]);
      sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a3.size(), x2.size());
      // tensor label: I1352
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x2, x1, x0, a3, c2, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x2, x1, x0, a3, c2, a1)]);
      sort_indices<1,4,0,2,3,5,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x2.size(), x1.size(), x0.size(), a3.size(), c2.size(), a1.size());
      dgemm_("T", "N", 1, ci0.size()*x1.size()*x0.size()*c2.size()*a1.size(), x2.size()*a3.size(),
             1.0, i0data_sorted, x2.size()*a3.size(), i1data_sorted, x2.size()*a3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,3,4,1,1,1,1>(odata_sorted, odata, ci0.size(), x1.size(), x0.size(), c2.size(), a1.size());
  out()->put_block(odata, ci0, x1, x0, c2, a1);
}

void Task1223::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  const Index a3 = b(4);
  const Index c2 = b(5);
  const Index a1 = b(6);
  // tensor label: I1352
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x2, x1, x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x2, x1, x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x2, x1, x0, a3, c2, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a1.size());
    // tensor label: I1353
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
    sort_indices<1,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), ci0.size()*x2.size()*x1.size()*x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }
  sort_indices<3,4,5,6,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), ci0.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x2, x1, x0, a3, c2, a1);
}

void Task1224::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1353
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,1,1,-1,4>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1225::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  const Index a3 = b(4);
  const Index c2 = b(5);
  const Index a1 = b(6);
  // tensor label: I1352
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x2, x1, x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x2, x1, x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x2, x1, x0, a3, c2, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());
    // tensor label: I1357
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x1, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x1, x2)]);
    sort_indices<1,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x1.size(), x2.size());
    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), ci0.size()*x0.size()*x1.size()*x2.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }
  sort_indices<3,6,5,4,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), ci0.size(), x0.size(), x1.size(), x2.size());
  out()->put_block(odata, ci0, x2, x1, x0, a3, c2, a1);
}

void Task1226::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  const Index x2 = b(4);
  // tensor label: I1357
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x1, x2);
  {
    // tensor label: Gamma410
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x1, x2);
    sort_indices<0,1,2,3,4,1,1,1,4>(i0data, odata, ci0.size(), x3.size(), x0.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, ci0, x3, x0, x1, x2);
}

void Task1227::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1305
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, x2)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), x2.size());
      // tensor label: I1934
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x1, x2)]);
      sort_indices<4,1,0,2,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x1.size(), x2.size());
      dgemm_("T", "N", a1.size()*c2.size(), ci0.size()*x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<2,4,3,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), ci0.size(), x0.size(), x1.size());
  out()->put_block(odata, ci0, x1, x0, c2, a1);
}

void Task1228::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  const Index x2 = b(4);
  // tensor label: I1934
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x1, x2);
  {
    // tensor label: Gamma410
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x1, x2);
    dscal_(ci0.size()*x3.size()*x0.size()*x1.size()*x2.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,4,1,1,-1,4>(i0data, odata, ci0.size(), x3.size(), x0.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, ci0, x3, x0, x1, x2);
}

void Task1229::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1305
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x3, x2)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x3.size(), x2.size());
      // tensor label: I1937
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
      sort_indices<2,1,0,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), ci0.size()*x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c2, a1);
}

void Task1230::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1937
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    dscal_(ci0.size()*x3.size()*x2.size()*x1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,4,1,1,1,4>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1231::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1305
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x3, x2)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x3.size(), x2.size());
      // tensor label: I2009
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
      sort_indices<2,1,0,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), ci0.size()*x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c2, a1);
}

void Task1232::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I2009
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1233::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1305
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x3, x2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x3, x2, a1)]);
      sort_indices<2,1,0,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), x3.size(), x2.size(), a1.size());
      // tensor label: I2012
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x3, x2, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x3, x2, x0)]);
      sort_indices<3,2,0,1,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x3.size(), x2.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), ci0.size()*x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c2, a1);
}

void Task1234::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x0 = b(4);
  // tensor label: I2012
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x3, x2, x0);
  {
    // tensor label: Gamma407
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x3, x2, x0);
    sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x1.size(), x3.size(), x2.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x3, x2, x0);
}

void Task1235::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1305
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, x2)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), x2.size());
      // tensor label: I2015
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x1, x2)]);
      sort_indices<4,1,0,2,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x1.size(), x2.size());
      dgemm_("T", "N", a1.size()*c2.size(), ci0.size()*x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<2,4,3,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), ci0.size(), x0.size(), x1.size());
  out()->put_block(odata, ci0, x1, x0, c2, a1);
}

void Task1236::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  const Index x2 = b(4);
  // tensor label: I2015
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x1, x2);
  {
    // tensor label: Gamma410
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x1, x2);
    sort_indices<0,1,2,3,4,1,1,1,2>(i0data, odata, ci0.size(), x3.size(), x0.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, ci0, x3, x0, x1, x2);
}

void Task1237::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1305
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, c2, a1)]);
      sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), c2.size(), a1.size());
      // tensor label: I2018
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
      sort_indices<2,1,0,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), ci0.size()*x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c2, a1);
}

void Task1238::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I2018
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1239::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1305
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a1), 0.0);
  // tensor label: h1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size());
  // tensor label: I2105
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
  sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
  dgemm_("T", "N", c2.size()*a1.size(), ci0.size()*x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a1.size());
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c2, a1);
}

void Task1240::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I2105
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1241::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& x0 : *range_[1]) {
        for (auto& x1 : *range_[1]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x0, x1)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x0.size(), x1.size());
          // tensor label: I1359
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0, c1, a2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0, c1, a2)]);
          sort_indices<3,4,2,1,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size(), c1.size(), a2.size());
          dgemm_("T", "N", 1, ci0.size(), x1.size()*x0.size()*c1.size()*a2.size(),
                 1.0, i0data_sorted, x1.size()*x0.size()*c1.size()*a2.size(), i1data_sorted, x1.size()*x0.size()*c1.size()*a2.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1242::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  const Index a2 = b(4);
  // tensor label: I1359
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c1, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), a2.size());
    // tensor label: I1360
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x2, x1, x0, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x2, x1, x0, c1)]);
    sort_indices<1,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x2.size(), x1.size(), x0.size(), c1.size());
    dgemm_("T", "N", a2.size(), ci0.size()*x1.size()*x0.size()*c1.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a2.size(), ci0.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, ci0, x1, x0, c1, a2);
}

void Task1243::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  const Index c1 = b(4);
  // tensor label: I1360
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x2, x1, x0, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x2, x1, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x2, x1, x0, c1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c1, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, c1, x3)]);
        sort_indices<3,1,0,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), c1.size(), x3.size());
        // tensor label: I1361
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x4, x2, x3, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x4, x2, x3, x1, x0)]);
        sort_indices<4,2,1,0,3,5,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
        dgemm_("T", "N", c1.size(), ci0.size()*x2.size()*x1.size()*x0.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, c1.size(), ci0.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x2, x1, x0, c1);
}

void Task1244::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x2 = b(3);
  const Index x3 = b(4);
  const Index x1 = b(5);
  const Index x0 = b(6);
  // tensor label: I1361
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x2, x3, x1, x0);
  {
    // tensor label: Gamma384
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x4, x2, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,6,1,1,1,4>(i0data, odata, ci0.size(), x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x5, x4, x2, x3, x1, x0);
}

void Task1245::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  const Index a2 = b(4);
  // tensor label: I1359
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, c3)]);
      sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x2.size(), c3.size());
      // tensor label: I1364
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x2, x1, x0, c3, a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x2, x1, x0, c3, a2, c1)]);
      sort_indices<4,1,0,2,3,5,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x2.size(), x1.size(), x0.size(), c3.size(), a2.size(), c1.size());
      dgemm_("T", "N", 1, ci0.size()*x1.size()*x0.size()*a2.size()*c1.size(), x2.size()*c3.size(),
             1.0, i0data_sorted, x2.size()*c3.size(), i1data_sorted, x2.size()*c3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,4,3,1,1,1,1>(odata_sorted, odata, ci0.size(), x1.size(), x0.size(), a2.size(), c1.size());
  out()->put_block(odata, ci0, x1, x0, c1, a2);
}

void Task1246::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  const Index c3 = b(4);
  const Index a2 = b(5);
  const Index c1 = b(6);
  // tensor label: I1364
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x2, x1, x0, c3, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x2, x1, x0, c3, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x2, x1, x0, c3, a2, c1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c1, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c1.size(), x3.size());
    // tensor label: I1365
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x2, x3, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x2, x3, x1, x0)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x2.size(), x3.size(), x1.size(), x0.size());
    dgemm_("T", "N", c3.size()*a2.size()*c1.size(), ci0.size()*x2.size()*x1.size()*x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c3.size()*a2.size()*c1.size());
  }
  sort_indices<3,4,5,6,0,1,2,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), ci0.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x2, x1, x0, c3, a2, c1);
}

void Task1247::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x2 = b(1);
  const Index x3 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1365
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x2, x3, x1, x0);
  {
    // tensor label: Gamma385
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x2, x3, x1, x0);
    sort_indices<0,1,2,3,4,1,1,-1,4>(i0data, odata, ci0.size(), x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x2, x3, x1, x0);
}

void Task1248::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  const Index c3 = b(4);
  const Index a2 = b(5);
  const Index c1 = b(6);
  // tensor label: I1364
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x2, x1, x0, c3, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x2, x1, x0, c3, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x2, x1, x0, c3, a2, c1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x3.size());
    // tensor label: I1369
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x2, x3, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x2, x3, x1, x0)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x2.size(), x3.size(), x1.size(), x0.size());
    dgemm_("T", "N", c1.size()*a2.size()*c3.size(), ci0.size()*x2.size()*x1.size()*x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<3,4,5,6,2,1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), ci0.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x2, x1, x0, c3, a2, c1);
}

void Task1249::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x2 = b(1);
  const Index x3 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1369
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x2, x3, x1, x0);
  {
    // tensor label: Gamma385
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x2, x3, x1, x0);
    sort_indices<0,1,2,3,4,1,1,1,2>(i0data, odata, ci0.size(), x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x2, x3, x1, x0);
}

