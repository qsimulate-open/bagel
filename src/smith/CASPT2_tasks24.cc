//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks24.cc
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


#include <src/smith/CASPT2_tasks24.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task1150::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  // tensor label: I1252
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x3);
  {
    // tensor label: Gamma392
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x3);
    sort_indices<0,1,2,1,1,-1,4>(i0data, odata, ci0.size(), x0.size(), x3.size());
  }
  out()->put_block(odata, ci0, x0, x3);
}

void Task1151::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);
  // tensor label: I1239
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c1, c3, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x3.size());
    // tensor label: I1255
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x3)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x3.size());
    dgemm_("T", "N", c1.size()*a2.size()*c3.size(), ci0.size()*x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<3,4,0,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, c1, c3, a2);
}

void Task1152::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  // tensor label: I1255
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x3);
  {
    // tensor label: Gamma392
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x3);
    sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x0.size(), x3.size());
  }
  out()->put_block(odata, ci0, x0, x3);
}

void Task1153::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);
  // tensor label: I1239
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c1, c3, a2), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, c4)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), c4.size());
    // tensor label: I1258
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, c4, a2, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, c4, a2, c3)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), c4.size(), a2.size(), c3.size());
    dgemm_("T", "N", c1.size(), ci0.size()*x0.size()*a2.size()*c3.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c1.size());
  }
  sort_indices<1,2,0,4,3,1,1,1,1>(odata_sorted, odata, c1.size(), ci0.size(), x0.size(), a2.size(), c3.size());
  out()->put_block(odata, ci0, x0, c1, c3, a2);
}

void Task1154::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c4 = b(2);
  const Index a2 = b(3);
  const Index c3 = b(4);
  // tensor label: I1258
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c4, a2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c4, a2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c4, a2, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a2, c3, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c3.size(), x1.size());
    // tensor label: I1259
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size());
    dgemm_("T", "N", c4.size()*a2.size()*c3.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c4.size()*a2.size()*c3.size());
  }
  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, c4.size(), a2.size(), c3.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, c4, a2, c3);
}

void Task1155::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  // tensor label: I1259
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1);
  {
    // tensor label: Gamma394
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
    sort_indices<0,1,2,1,1,-1,2>(i0data, odata, ci0.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x0, x1);
}

void Task1156::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c4 = b(2);
  const Index a2 = b(3);
  const Index c3 = b(4);
  // tensor label: I1258
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c4, a2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c4, a2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c4, a2, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c4, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c4, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c4.size(), x1.size());
    // tensor label: I1263
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size());
    dgemm_("T", "N", c3.size()*a2.size()*c4.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c3.size()*a2.size()*c4.size());
  }
  sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c4.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, c4, a2, c3);
}

void Task1157::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  // tensor label: I1263
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1);
  {
    // tensor label: Gamma394
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
    sort_indices<0,1,2,1,1,1,4>(i0data, odata, ci0.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x0, x1);
}

void Task1158::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);
  // tensor label: I1239
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c1, c3, a2), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, c4)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), c4.size());
    // tensor label: I1266
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, c4, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, c4, a2, c1)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), c4.size(), a2.size(), c1.size());
    dgemm_("T", "N", c3.size(), ci0.size()*x0.size()*a2.size()*c1.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c3.size());
  }
  sort_indices<1,2,4,0,3,1,1,1,1>(odata_sorted, odata, c3.size(), ci0.size(), x0.size(), a2.size(), c1.size());
  out()->put_block(odata, ci0, x0, c1, c3, a2);
}

void Task1159::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c4 = b(2);
  const Index a2 = b(3);
  const Index c1 = b(4);
  // tensor label: I1266
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c4, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c4, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c4, a2, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c1, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a2, c1, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c1.size(), x1.size());
    // tensor label: I1267
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size());
    dgemm_("T", "N", c4.size()*a2.size()*c1.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c4.size()*a2.size()*c1.size());
  }
  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, c4.size(), a2.size(), c1.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, c4, a2, c1);
}

void Task1160::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  // tensor label: I1267
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1);
  {
    // tensor label: Gamma394
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
    sort_indices<0,1,2,1,1,1,4>(i0data, odata, ci0.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x0, x1);
}

void Task1161::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c4 = b(2);
  const Index a2 = b(3);
  const Index c1 = b(4);
  // tensor label: I1266
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c4, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c4, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c4, a2, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c4, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c4, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c4.size(), x1.size());
    // tensor label: I1275
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size());
    dgemm_("T", "N", c1.size()*a2.size()*c4.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c4.size());
  }
  sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c4.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, c4, a2, c1);
}

void Task1162::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  // tensor label: I1275
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1);
  {
    // tensor label: Gamma394
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
    sort_indices<0,1,2,1,1,-1,2>(i0data, odata, ci0.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x0, x1);
}

void Task1163::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);
  // tensor label: I1239
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c1, c3, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a2.size());
    // tensor label: I1270
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, c3, a4, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, c3, a4, c1)]);
    sort_indices<3,0,1,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), c3.size(), a4.size(), c1.size());
    dgemm_("T", "N", a2.size(), ci0.size()*x0.size()*c3.size()*c1.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<1,2,4,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), ci0.size(), x0.size(), c3.size(), c1.size());
  out()->put_block(odata, ci0, x0, c1, c3, a2);
}

void Task1164::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c3 = b(2);
  const Index a4 = b(3);
  const Index c1 = b(4);
  // tensor label: I1270
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c3, a4, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c3, a4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c3, a4, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, c1, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4, c1, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), c1.size(), x1.size());
    // tensor label: I1271
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size());
    dgemm_("T", "N", c3.size()*a4.size()*c1.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c3.size()*a4.size()*c1.size());
  }
  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, c3.size(), a4.size(), c1.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, c3, a4, c1);
}

void Task1165::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  // tensor label: I1271
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1);
  {
    // tensor label: Gamma394
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
    sort_indices<0,1,2,1,1,-1,4>(i0data, odata, ci0.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x0, x1);
}

void Task1166::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c3 = b(2);
  const Index a4 = b(3);
  const Index c1 = b(4);
  // tensor label: I1270
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c3, a4, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c3, a4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c3, a4, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x1.size());
    // tensor label: I1279
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size());
    dgemm_("T", "N", c1.size()*a4.size()*c3.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size()*a4.size()*c3.size());
  }
  sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, c3, a4, c1);
}

void Task1167::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  // tensor label: I1279
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1);
  {
    // tensor label: Gamma394
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
    sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x0, x1);
}

void Task1168::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);
  // tensor label: I1239
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c1, c3, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), x1.size());
    // tensor label: I1282
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0, a2, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0, a2, c3)]);
    sort_indices<1,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size(), a2.size(), c3.size());
    dgemm_("T", "N", c1.size(), ci0.size()*x0.size()*a2.size()*c3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size());
  }
  sort_indices<1,2,0,4,3,1,1,1,1>(odata_sorted, odata, c1.size(), ci0.size(), x0.size(), a2.size(), c3.size());
  out()->put_block(odata, ci0, x0, c1, c3, a2);
}

void Task1169::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  const Index c3 = b(4);
  // tensor label: I1282
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, a2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, a2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, a2, c3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c3, x2)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c3.size(), x2.size());
      // tensor label: I1283
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x1, x0, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x1, x0, x2)]);
      sort_indices<4,1,0,2,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x1.size(), x0.size(), x2.size());
      dgemm_("T", "N", a2.size()*c3.size(), ci0.size()*x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a2.size()*c3.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a2.size(), c3.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, a2, c3);
}

void Task1170::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  const Index x2 = b(4);
  // tensor label: I1283
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x1, x0, x2);
  {
    // tensor label: Gamma400
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x1, x0, x2);
    sort_indices<0,1,2,3,4,1,1,1,4>(i0data, odata, ci0.size(), x3.size(), x1.size(), x0.size(), x2.size());
  }
  out()->put_block(odata, ci0, x3, x1, x0, x2);
}

void Task1171::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  const Index c3 = b(4);
  // tensor label: I1282
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, a2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, a2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, a2, c3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, x3, x2)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), x3.size(), x2.size());
      // tensor label: I1291
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x0, x1)]);
      sort_indices<2,1,0,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x0.size(), x1.size());
      dgemm_("T", "N", c3.size()*a2.size(), ci0.size()*x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c3.size()*a2.size());
    }
  }
  sort_indices<2,4,3,1,0,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), ci0.size(), x0.size(), x1.size());
  out()->put_block(odata, ci0, x1, x0, a2, c3);
}

void Task1172::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x0 = b(3);
  const Index x1 = b(4);
  // tensor label: I1291
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x0, x1);
  {
    // tensor label: Gamma390
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x0, x1);
    sort_indices<0,1,2,3,4,1,1,-1,4>(i0data, odata, ci0.size(), x3.size(), x2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x2, x0, x1);
}

void Task1173::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);
  // tensor label: I1239
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c1, c3, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), x1.size());
    // tensor label: I1286
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1, a2, c1)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size(), a2.size(), c1.size());
    dgemm_("T", "N", c3.size(), ci0.size()*x0.size()*a2.size()*c1.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c3.size());
  }
  sort_indices<1,2,4,0,3,1,1,1,1>(odata_sorted, odata, c3.size(), ci0.size(), x0.size(), a2.size(), c1.size());
  out()->put_block(odata, ci0, x0, c1, c3, a2);
}

void Task1174::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  const Index c1 = b(4);
  // tensor label: I1286
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a2, c1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c1, x2)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c1.size(), x2.size());
      // tensor label: I1287
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x0, x1)]);
      sort_indices<2,1,0,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x0.size(), x1.size());
      dgemm_("T", "N", a2.size()*c1.size(), ci0.size()*x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), ci0.size(), x0.size(), x1.size());
  out()->put_block(odata, ci0, x0, x1, a2, c1);
}

void Task1175::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x0 = b(3);
  const Index x1 = b(4);
  // tensor label: I1287
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x0, x1);
  {
    // tensor label: Gamma390
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x0, x1);
    sort_indices<0,1,2,3,4,1,1,-1,4>(i0data, odata, ci0.size(), x3.size(), x2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x2, x0, x1);
}

void Task1176::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  const Index c1 = b(4);
  // tensor label: I1286
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a2, c1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x3, x2)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x3.size(), x2.size());
      // tensor label: I1295
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x0, x1)]);
      sort_indices<2,1,0,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x0.size(), x1.size());
      dgemm_("T", "N", c1.size()*a2.size(), ci0.size()*x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<2,3,4,1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), ci0.size(), x0.size(), x1.size());
  out()->put_block(odata, ci0, x0, x1, a2, c1);
}

void Task1177::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x0 = b(3);
  const Index x1 = b(4);
  // tensor label: I1295
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x0, x1);
  {
    // tensor label: Gamma390
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x0, x1);
    sort_indices<0,1,2,3,4,1,1,1,2>(i0data, odata, ci0.size(), x3.size(), x2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x2, x0, x1);
}

void Task1178::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);
  // tensor label: I1239
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c1, c3, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, x1)]);
      sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a4.size(), x1.size());
      // tensor label: I1298
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1, c1, a4, c3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1, c1, a4, c3, a2)]);
      sort_indices<2,4,0,1,3,5,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size(), c1.size(), a4.size(), c3.size(), a2.size());
      dgemm_("T", "N", 1, ci0.size()*x0.size()*c1.size()*c3.size()*a2.size(), x1.size()*a4.size(),
             1.0, i0data_sorted, x1.size()*a4.size(), i1data_sorted, x1.size()*a4.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,3,4,1,1,1,1>(odata_sorted, odata, ci0.size(), x0.size(), c1.size(), c3.size(), a2.size());
  out()->put_block(odata, ci0, x0, c1, c3, a2);
}

void Task1179::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index c1 = b(3);
  const Index a4 = b(4);
  const Index c3 = b(5);
  const Index a2 = b(6);
  // tensor label: I1298
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, c1, a4, c3, a2), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
  // tensor label: I1299
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1)]);
  sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size());
  dgemm_("T", "N", c1.size()*a4.size()*c3.size()*a2.size(), ci0.size()*x0.size()*x1.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a4.size()*c3.size()*a2.size());
  sort_indices<4,5,6,0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size(), ci0.size(), x0.size(), x1.size());
  out()->put_block(odata, ci0, x0, x1, c1, a4, c3, a2);
}

void Task1180::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  // tensor label: I1299
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1);
  {
    // tensor label: Gamma394
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
    sort_indices<0,1,2,1,1,-1,2>(i0data, odata, ci0.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x0, x1);
}

void Task1181::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index c1 = b(3);
  const Index a4 = b(4);
  const Index c3 = b(5);
  const Index a2 = b(6);
  // tensor label: I1298
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, c1, a4, c3, a2), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
  // tensor label: I1303
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1)]);
  sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size());
  dgemm_("T", "N", c1.size()*a2.size()*c3.size()*a4.size(), ci0.size()*x0.size()*x1.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a2.size()*c3.size()*a4.size());
  sort_indices<4,5,6,0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size(), ci0.size(), x0.size(), x1.size());
  out()->put_block(odata, ci0, x0, x1, c1, a4, c3, a2);
}

void Task1182::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  // tensor label: I1303
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1);
  {
    // tensor label: Gamma394
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
    sort_indices<0,1,2,1,1,1,1>(i0data, odata, ci0.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x0, x1);
}

void Task1183::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);
  // tensor label: I1239
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c1, c3, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c1, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c1.size(), x1.size());
    // tensor label: I1928
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size());
    dgemm_("T", "N", c3.size()*a2.size()*c1.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c3.size()*a2.size()*c1.size());
  }
  sort_indices<3,4,2,0,1,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, c1, c3, a2);
}

void Task1184::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  // tensor label: I1928
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1);
  {
    // tensor label: Gamma394
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
    dscal_(ci0.size()*x0.size()*x1.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,1,1,1,4>(i0data, odata, ci0.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x0, x1);
}

void Task1185::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);
  // tensor label: I1239
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c1, c3, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x1.size());
    // tensor label: I1931
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size());
    dgemm_("T", "N", c1.size()*a2.size()*c3.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<3,4,0,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, c1, c3, a2);
}

void Task1186::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  // tensor label: I1931
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1);
  {
    // tensor label: Gamma394
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
    dscal_(ci0.size()*x0.size()*x1.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,1,1,-1,2>(i0data, odata, ci0.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x0, x1);
}

void Task1187::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);
  // tensor label: I1239
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c1, c3, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x1, c1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x1, c1, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), x1.size(), c1.size(), a2.size());
    // tensor label: I2003
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size());
    dgemm_("T", "N", c3.size()*c1.size()*a2.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c3.size()*c1.size()*a2.size());
  }
  sort_indices<3,4,1,0,2,1,1,1,1>(odata_sorted, odata, c3.size(), c1.size(), a2.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, c1, c3, a2);
}

void Task1188::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  // tensor label: I2003
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1);
  {
    // tensor label: Gamma394
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
    sort_indices<0,1,2,1,1,2,1>(i0data, odata, ci0.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x0, x1);
}

void Task1189::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);
  // tensor label: I1239
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c1, c3, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x1, c3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x1, c3, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x1.size(), c3.size(), a2.size());
    // tensor label: I2006
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size());
    dgemm_("T", "N", c1.size()*c3.size()*a2.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size()*c3.size()*a2.size());
  }
  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), a2.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, c1, c3, a2);
}

void Task1190::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  // tensor label: I2006
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1);
  {
    // tensor label: Gamma394
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
    sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x0, x1);
}

void Task1191::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& c2 : *range_[0]) {
        for (auto& x1 : *range_[1]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, x1)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), x1.size());
          // tensor label: I1305
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0, c2, a1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0, c2, a1)]);
          sort_indices<2,4,3,1,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size(), c2.size(), a1.size());
          dgemm_("T", "N", 1, ci0.size(), x1.size()*x0.size()*c2.size()*a1.size(),
                 1.0, i0data_sorted, x1.size()*x0.size()*c2.size()*a1.size(), i1data_sorted, x1.size()*x0.size()*c2.size()*a1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1192::Task_local::compute() {
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
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), a1.size());
    // tensor label: I1306
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x2, x0, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x2, x0, c2)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x2.size(), x0.size(), c2.size());
    dgemm_("T", "N", a1.size(), ci0.size()*x1.size()*x0.size()*c2.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x1.size(), x0.size(), c2.size());
  out()->put_block(odata, ci0, x1, x0, c2, a1);
}

void Task1193::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index x0 = b(3);
  const Index c2 = b(4);
  // tensor label: I1306
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x2, x0, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x2, x0, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x2, x0, c2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c2, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, c2, x3)]);
        sort_indices<3,1,0,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), c2.size(), x3.size());
        // tensor label: I1307
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x4, x1, x3, x2, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x4, x1, x3, x2, x0)]);
        sort_indices<4,2,1,0,3,5,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x4.size(), x1.size(), x3.size(), x2.size(), x0.size());
        dgemm_("T", "N", c2.size(), ci0.size()*x1.size()*x2.size()*x0.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, c2.size());
      }
    }
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, c2.size(), ci0.size(), x1.size(), x2.size(), x0.size());
  out()->put_block(odata, ci0, x1, x2, x0, c2);
}

void Task1194::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x1 = b(3);
  const Index x3 = b(4);
  const Index x2 = b(5);
  const Index x0 = b(6);
  // tensor label: I1307
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x1, x3, x2, x0);
  {
    // tensor label: Gamma406
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x4, x1, x3, x2, x0);
    sort_indices<0,1,2,3,4,5,6,1,1,-1,4>(i0data, odata, ci0.size(), x5.size(), x4.size(), x1.size(), x3.size(), x2.size(), x0.size());
  }
  out()->put_block(odata, ci0, x5, x4, x1, x3, x2, x0);
}

void Task1195::Task_local::compute() {
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
    for (auto& x2 : *range_[1]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, c3)]);
      sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x2.size(), c3.size());
      // tensor label: I1310
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x2, x0, c3, a1, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x2, x0, c3, a1, c2)]);
      sort_indices<4,2,0,1,3,5,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x2.size(), x0.size(), c3.size(), a1.size(), c2.size());
      dgemm_("T", "N", 1, ci0.size()*x1.size()*x0.size()*a1.size()*c2.size(), x2.size()*c3.size(),
             1.0, i0data_sorted, x2.size()*c3.size(), i1data_sorted, x2.size()*c3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,4,3,1,1,1,1>(odata_sorted, odata, ci0.size(), x1.size(), x0.size(), a1.size(), c2.size());
  out()->put_block(odata, ci0, x1, x0, c2, a1);
}

void Task1196::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index x0 = b(3);
  const Index c3 = b(4);
  const Index a1 = b(5);
  const Index c2 = b(6);
  // tensor label: I1310
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x2, x0, c3, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x2, x0, c3, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x2, x0, c3, a1, c2), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a1, c2, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a1, c2, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), a1.size(), c2.size(), x3.size());
    // tensor label: I1311
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x3, x2, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x3, x2, x0)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x3.size(), x2.size(), x0.size());
    dgemm_("T", "N", c3.size()*a1.size()*c2.size(), ci0.size()*x1.size()*x2.size()*x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c3.size()*a1.size()*c2.size());
  }
  sort_indices<3,4,5,6,0,1,2,1,1,1,1>(odata_sorted, odata, c3.size(), a1.size(), c2.size(), ci0.size(), x1.size(), x2.size(), x0.size());
  out()->put_block(odata, ci0, x1, x2, x0, c3, a1, c2);
}

void Task1197::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x0 = b(4);
  // tensor label: I1311
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x3, x2, x0);
  {
    // tensor label: Gamma407
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x3, x2, x0);
    sort_indices<0,1,2,3,4,1,1,1,4>(i0data, odata, ci0.size(), x1.size(), x3.size(), x2.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x3, x2, x0);
}

void Task1198::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index x0 = b(3);
  const Index c3 = b(4);
  const Index a1 = b(5);
  const Index c2 = b(6);
  // tensor label: I1310
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x2, x0, c3, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x2, x0, c3, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x2, x0, c3, a1, c2), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c3, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c3, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c3.size(), x3.size());
    // tensor label: I1315
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x2, x3, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x2, x3, x1, x0)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x2.size(), x3.size(), x1.size(), x0.size());
    dgemm_("T", "N", c2.size()*a1.size()*c3.size(), ci0.size()*x2.size()*x1.size()*x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c2.size()*a1.size()*c3.size());
  }
  sort_indices<3,5,4,6,2,1,0,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), c3.size(), ci0.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x2, x0, c3, a1, c2);
}

void Task1199::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x2 = b(1);
  const Index x3 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1315
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x2, x3, x1, x0);
  {
    // tensor label: Gamma385
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x2, x3, x1, x0);
    sort_indices<0,1,2,3,4,1,1,-1,4>(i0data, odata, ci0.size(), x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x2, x3, x1, x0);
}

