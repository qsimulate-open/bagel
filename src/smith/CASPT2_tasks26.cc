//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks26.cc
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


#include <src/smith/CASPT2_tasks26.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task1250::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  const Index a2 = b(4);
  // tensor label: I1359
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c1, a2), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, c1, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, c1, x4)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), c1.size(), x4.size());
      // tensor label: I1372
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x4, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x4, x1, x0)]);
      sort_indices<2,1,0,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x4.size(), x1.size(), x0.size());
      dgemm_("T", "N", a2.size()*c1.size(), ci0.size()*x1.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<2,3,4,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c1, a2);
}

void Task1251::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1372
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x1, x0);
  {
    // tensor label: Gamma412
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x4, x1, x0);
    sort_indices<0,1,2,3,4,1,1,-1,4>(i0data, odata, ci0.size(), x5.size(), x4.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x5, x4, x1, x0);
}

void Task1252::Task_local::compute() {
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
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, c3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), c3.size());
    // tensor label: I1375
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0, a2, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0, a2, c3)]);
    sort_indices<4,0,1,2,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size(), a2.size(), c3.size());
    dgemm_("T", "N", c1.size(), ci0.size()*x1.size()*x0.size()*a2.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, c1.size());
  }
  sort_indices<1,2,3,0,4,1,1,1,1>(odata_sorted, odata, c1.size(), ci0.size(), x1.size(), x0.size(), a2.size());
  out()->put_block(odata, ci0, x1, x0, c1, a2);
}

void Task1253::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  const Index c3 = b(4);
  // tensor label: I1375
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, a2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, a2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, a2, c3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c3, x2)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c3.size(), x2.size());
      // tensor label: I1376
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
      sort_indices<2,1,0,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", a2.size()*c3.size(), ci0.size()*x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a2.size()*c3.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a2.size(), c3.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, a2, c3);
}

void Task1254::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1376
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,1,1,1,4>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1255::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  const Index c3 = b(4);
  // tensor label: I1375
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, a2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, a2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, a2, c3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, x3, x2)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), x3.size(), x2.size());
      // tensor label: I1387
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
      sort_indices<2,1,0,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c3.size()*a2.size(), ci0.size()*x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c3.size()*a2.size());
    }
  }
  sort_indices<2,3,4,1,0,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, a2, c3);
}

void Task1256::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1387
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1257::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  const Index a2 = b(4);
  // tensor label: I1359
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c1, a2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), a2.size());
    // tensor label: I1379
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0, a3, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0, a3, c1)]);
    sort_indices<3,0,1,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size(), a3.size(), c1.size());
    dgemm_("T", "N", a2.size(), ci0.size()*x1.size()*x0.size()*c1.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a2.size(), ci0.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, ci0, x1, x0, c1, a2);
}

void Task1258::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a3 = b(3);
  const Index c1 = b(4);
  // tensor label: I1379
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, a3, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, a3, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, a3, c1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c1, x2)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c1.size(), x2.size());
      // tensor label: I1380
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
      sort_indices<2,1,0,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", a3.size()*c1.size(), ci0.size()*x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a3.size()*c1.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a3.size(), c1.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, a3, c1);
}

void Task1259::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1380
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,1,1,-1,4>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1260::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a3 = b(3);
  const Index c1 = b(4);
  // tensor label: I1379
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, a3, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, a3, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, a3, c1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, x3, x2)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), x3.size(), x2.size());
      // tensor label: I1391
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
      sort_indices<2,1,0,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a3.size(), ci0.size()*x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c1.size()*a3.size());
    }
  }
  sort_indices<2,3,4,1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a3.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, a3, c1);
}

void Task1261::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1391
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,1,1,1,2>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1262::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  const Index a2 = b(4);
  // tensor label: I1359
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c1, a2), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x5, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x5, x4)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x5.size(), x4.size());
      // tensor label: I1383
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x4, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x4, x1, x0)]);
      sort_indices<2,1,0,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x4.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), ci0.size()*x1.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c1, a2);
}

void Task1263::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1383
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x1, x0);
  {
    // tensor label: Gamma412
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x4, x1, x0);
    sort_indices<0,1,2,3,4,1,1,1,2>(i0data, odata, ci0.size(), x5.size(), x4.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x5, x4, x1, x0);
}

void Task1264::Task_local::compute() {
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
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), x2.size());
    // tensor label: I1394
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x2, x1, x0, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x2, x1, x0, a2)]);
    sort_indices<1,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x2.size(), x1.size(), x0.size(), a2.size());
    dgemm_("T", "N", c1.size(), ci0.size()*x1.size()*x0.size()*a2.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, c1.size());
  }
  sort_indices<1,2,3,0,4,1,1,1,1>(odata_sorted, odata, c1.size(), ci0.size(), x1.size(), x0.size(), a2.size());
  out()->put_block(odata, ci0, x1, x0, c1, a2);
}

void Task1265::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  const Index a2 = b(4);
  // tensor label: I1394
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x2, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x2, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x2, x1, x0, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, x4, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, x4, x3)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), x4.size(), x3.size());
        // tensor label: I1395
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x2, x4, x3, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x2, x4, x3, x1, x0)]);
        sort_indices<4,3,1,0,2,5,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
        dgemm_("T", "N", a2.size(), ci0.size()*x2.size()*x1.size()*x0.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a2.size(), ci0.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x2, x1, x0, a2);
}

void Task1266::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x2 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x1 = b(5);
  const Index x0 = b(6);
  // tensor label: I1395
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x2, x4, x3, x1, x0);
  {
    // tensor label: Gamma429
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x2, x4, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,6,1,1,-1,4>(i0data, odata, ci0.size(), x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x5, x2, x4, x3, x1, x0);
}

void Task1267::Task_local::compute() {
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
    for (auto& a4 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, c3)]);
      sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a4.size(), c3.size());
      // tensor label: I1398
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0, c1, a4, c3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0, c1, a4, c3, a2)]);
      sort_indices<5,4,0,1,2,3,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size(), c1.size(), a4.size(), c3.size(), a2.size());
      dgemm_("T", "N", 1, ci0.size()*x1.size()*x0.size()*c1.size()*a2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,3,4,1,1,1,1>(odata_sorted, odata, ci0.size(), x1.size(), x0.size(), c1.size(), a2.size());
  out()->put_block(odata, ci0, x1, x0, c1, a2);
}

void Task1268::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  const Index a4 = b(4);
  const Index c3 = b(5);
  const Index a2 = b(6);
  // tensor label: I1398
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c1, a4, c3, a2), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
  // tensor label: I1399
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
  sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
  dgemm_("T", "N", c1.size()*a4.size()*c3.size()*a2.size(), ci0.size()*x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a4.size()*c3.size()*a2.size());
  sort_indices<4,5,6,0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c1, a4, c3, a2);
}

void Task1269::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1399
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1270::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  const Index a4 = b(4);
  const Index c3 = b(5);
  const Index a2 = b(6);
  // tensor label: I1398
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c1, a4, c3, a2), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
  // tensor label: I1403
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
  sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
  dgemm_("T", "N", c1.size()*a2.size()*c3.size()*a4.size(), ci0.size()*x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a2.size()*c3.size()*a4.size());
  sort_indices<4,5,6,0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c1, a4, c3, a2);
}

void Task1271::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1403
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,2,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1272::Task_local::compute() {
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
    for (auto& a3 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, x2)]);
      sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a3.size(), x2.size());
      // tensor label: I1406
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x2, x1, x0, a3, c1, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x2, x1, x0, a3, c1, a2)]);
      sort_indices<1,4,0,2,3,5,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x2.size(), x1.size(), x0.size(), a3.size(), c1.size(), a2.size());
      dgemm_("T", "N", 1, ci0.size()*x1.size()*x0.size()*c1.size()*a2.size(), x2.size()*a3.size(),
             1.0, i0data_sorted, x2.size()*a3.size(), i1data_sorted, x2.size()*a3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,3,4,1,1,1,1>(odata_sorted, odata, ci0.size(), x1.size(), x0.size(), c1.size(), a2.size());
  out()->put_block(odata, ci0, x1, x0, c1, a2);
}

void Task1273::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  const Index a3 = b(4);
  const Index c1 = b(5);
  const Index a2 = b(6);
  // tensor label: I1406
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x2, x1, x0, a3, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x2, x1, x0, a3, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x2, x1, x0, a3, c1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c1, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c1.size(), a2.size());
    // tensor label: I1407
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
    sort_indices<1,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
    dgemm_("T", "N", a3.size()*c1.size()*a2.size(), ci0.size()*x2.size()*x1.size()*x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a3.size()*c1.size()*a2.size());
  }
  sort_indices<3,4,5,6,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c1.size(), a2.size(), ci0.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x2, x1, x0, a3, c1, a2);
}

void Task1274::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1407
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,1,1,1,2>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1275::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  const Index a3 = b(4);
  const Index c1 = b(5);
  const Index a2 = b(6);
  // tensor label: I1406
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x2, x1, x0, a3, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x2, x1, x0, a3, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x2, x1, x0, a3, c1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c1, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c1.size(), a3.size());
    // tensor label: I1411
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
    sort_indices<1,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
    dgemm_("T", "N", a2.size()*c1.size()*a3.size(), ci0.size()*x2.size()*x1.size()*x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a2.size()*c1.size()*a3.size());
  }
  sort_indices<3,4,5,6,2,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), a3.size(), ci0.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x2, x1, x0, a3, c1, a2);
}

void Task1276::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1411
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,1,1,-1,4>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1277::Task_local::compute() {
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
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c1, x2)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c1.size(), x2.size());
      // tensor label: I1940
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
      sort_indices<2,1,0,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", a2.size()*c1.size(), ci0.size()*x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<2,3,4,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c1, a2);
}

void Task1278::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1940
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    dscal_(ci0.size()*x3.size()*x2.size()*x1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,4,1,1,1,4>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1279::Task_local::compute() {
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
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x3, x2)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x3.size(), x2.size());
      // tensor label: I1943
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
      sort_indices<2,1,0,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), ci0.size()*x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c1, a2);
}

void Task1280::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1943
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    dscal_(ci0.size()*x3.size()*x2.size()*x1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1281::Task_local::compute() {
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
    for (auto& x3 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x3, x2)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x3.size(), x2.size());
      // tensor label: I2021
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
      sort_indices<2,1,0,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), ci0.size()*x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c1, a2);
}

void Task1282::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I2021
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1283::Task_local::compute() {
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
    for (auto& x3 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x3, x2, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x3, x2, a2)]);
      sort_indices<2,1,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x3.size(), x2.size(), a2.size());
      // tensor label: I2024
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x2, x3, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x2, x3, x1, x0)]);
      sort_indices<1,2,0,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x2.size(), x3.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), ci0.size()*x1.size()*x0.size(), x2.size()*x3.size(),
             1.0, i0data_sorted, x2.size()*x3.size(), i1data_sorted, x2.size()*x3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c1, a2);
}

void Task1284::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x2 = b(1);
  const Index x3 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I2024
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x2, x3, x1, x0);
  {
    // tensor label: Gamma385
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x2, x3, x1, x0);
    sort_indices<0,1,2,3,4,1,1,1,2>(i0data, odata, ci0.size(), x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x2, x3, x1, x0);
}

void Task1285::Task_local::compute() {
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
    for (auto& x3 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c1, x2)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c1.size(), x2.size());
      // tensor label: I2027
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
      sort_indices<2,1,0,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", a2.size()*c1.size(), ci0.size()*x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<2,3,4,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c1, a2);
}

void Task1286::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I2027
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1287::Task_local::compute() {
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
    for (auto& x3 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, c1, a2)]);
      sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), c1.size(), a2.size());
      // tensor label: I2030
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
      sort_indices<2,1,0,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), ci0.size()*x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c1, a2);
}

void Task1288::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I2030
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1289::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  const Index a2 = b(4);
  // tensor label: I1359
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c1, a2), 0.0);
  // tensor label: h1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size());
  // tensor label: I2108
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
  sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
  dgemm_("T", "N", c1.size()*a2.size(), ci0.size()*x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a2.size());
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c1, a2);
}

void Task1290::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I2108
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,2,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1291::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x1 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
          // tensor label: I1413
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x2, x1, a1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x2, x1, a1)]);
          sort_indices<1,4,3,2,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x2.size(), x1.size(), a1.size());
          dgemm_("T", "N", 1, ci0.size(), x0.size()*x2.size()*x1.size()*a1.size(),
                 1.0, i0data_sorted, x0.size()*x2.size()*x1.size()*a1.size(), i1data_sorted, x0.size()*x2.size()*x1.size()*a1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1292::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a1 = b(4);
  // tensor label: I1413
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x2, x1, a1), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, c2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, c2)]);
      sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), c2.size());
      // tensor label: I1414
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x3, x2, x1, a1, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x3, x2, x1, a1, c2)]);
      sort_indices<6,2,0,1,3,4,5,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x3.size(), x2.size(), x1.size(), a1.size(), c2.size());
      dgemm_("T", "N", 1, ci0.size()*x0.size()*x2.size()*x1.size()*a1.size(), x3.size()*c2.size(),
             1.0, i0data_sorted, x3.size()*c2.size(), i1data_sorted, x3.size()*c2.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,3,4,1,1,1,1>(odata_sorted, odata, ci0.size(), x0.size(), x2.size(), x1.size(), a1.size());
  out()->put_block(odata, ci0, x0, x2, x1, a1);
}

void Task1293::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  const Index a1 = b(5);
  const Index c2 = b(6);
  // tensor label: I1414
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x3, x2, x1, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x3, x2, x1, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x3, x2, x1, a1, c2), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c2, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c2, x4)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c2.size(), x4.size());
      // tensor label: I1415
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x0, x3, x4, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x0, x3, x4, x2, x1)]);
      sort_indices<4,1,0,2,3,5,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x0.size(), x3.size(), x4.size(), x2.size(), x1.size());
      dgemm_("T", "N", a1.size()*c2.size(), ci0.size()*x0.size()*x3.size()*x2.size()*x1.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<2,3,4,5,6,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), ci0.size(), x0.size(), x3.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x0, x3, x2, x1, a1, c2);
}

void Task1294::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  const Index x4 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);
  // tensor label: I1415
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x3, x4, x2, x1);
  {
    // tensor label: Gamma434
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x3, x4, x2, x1);
    sort_indices<0,1,2,3,4,5,6,1,1,1,4>(i0data, odata, ci0.size(), x5.size(), x0.size(), x3.size(), x4.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x0, x3, x4, x2, x1);
}

void Task1295::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  const Index a1 = b(5);
  const Index c2 = b(6);
  // tensor label: I1414
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x3, x2, x1, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x3, x2, x1, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x3, x2, x1, a1, c2), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x5, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x5, x4)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x5.size(), x4.size());
      // tensor label: I1419
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x4, x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x4, x3, x0, x2, x1)]);
      sort_indices<2,1,0,3,4,5,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
      dgemm_("T", "N", c2.size()*a1.size(), ci0.size()*x3.size()*x0.size()*x2.size()*x1.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,4,3,5,6,1,0,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x0, x3, x2, x1, a1, c2);
}

void Task1296::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x0 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);
  // tensor label: I1419
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x3, x0, x2, x1);
  {
    // tensor label: Gamma435
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x4, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,5,6,1,1,-1,4>(i0data, odata, ci0.size(), x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x4, x3, x0, x2, x1);
}

void Task1297::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a1 = b(4);
  // tensor label: I1413
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x2, x1, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x7 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a1, x6, x5);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a1, x6, x5)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x7.size(), a1.size(), x6.size(), x5.size());
        // tensor label: I1422
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x7, x0, x6, x5, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x7, x0, x6, x5, x2, x1)]);
        sort_indices<4,3,1,0,2,5,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x7.size(), x0.size(), x6.size(), x5.size(), x2.size(), x1.size());
        dgemm_("T", "N", a1.size(), ci0.size()*x0.size()*x2.size()*x1.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x0, x2, x1, a1);
}

void Task1298::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x7 = b(1);
  const Index x0 = b(2);
  const Index x6 = b(3);
  const Index x5 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);
  // tensor label: I1422
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x7, x0, x6, x5, x2, x1);
  {
    // tensor label: Gamma436
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x7, x0, x6, x5, x2, x1);
    sort_indices<0,1,2,3,4,5,6,1,1,1,4>(i0data, odata, ci0.size(), x7.size(), x0.size(), x6.size(), x5.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x7, x0, x6, x5, x2, x1);
}

void Task1299::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a1 = b(4);
  // tensor label: I1413
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x2, x1, a1), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a2, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a2.size(), a1.size());
    // tensor label: I1425
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x2, x1, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x2, x1, a2)]);
    sort_indices<4,0,1,2,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x2.size(), x1.size(), a2.size());
    dgemm_("T", "N", a1.size(), ci0.size()*x0.size()*x2.size()*x1.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x0, x2, x1, a1);
}

