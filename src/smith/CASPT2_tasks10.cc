//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks10.cc
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


#include <src/smith/CASPT2_tasks10.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task450::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c1 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I391
  std::unique_ptr<double[]> odata = out()->move_block(x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c1, c3, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x3.size());
    // tensor label: I407
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x3.size());
    dgemm_("T", "N", c1.size()*a2.size()*c3.size(), x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<3,0,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), x0.size());
  out()->put_block(odata, x0, c1, c3, a2);
}

void Task451::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x3 = b(1);
  // tensor label: I407
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3);
  {
    // tensor label: Gamma14
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3);
    sort_indices<0,1,1,1,1,2>(i0data, odata, x0.size(), x3.size());
  }
  out()->put_block(odata, x0, x3);
}

void Task452::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c1 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I391
  std::unique_ptr<double[]> odata = out()->move_block(x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c1, c3, a2), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, c4)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), c4.size());
    // tensor label: I410
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c4, a2, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c4, a2, c3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c4.size(), a2.size(), c3.size());
    dgemm_("T", "N", c1.size(), x0.size()*a2.size()*c3.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c1.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), x0.size(), a2.size(), c3.size());
  out()->put_block(odata, x0, c1, c3, a2);
}

void Task453::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c4 = b(1);
  const Index a2 = b(2);
  const Index c3 = b(3);
  // tensor label: I410
  std::unique_ptr<double[]> odata = out()->move_block(x0, c4, a2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c4, a2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c4, a2, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a2, c3, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c3.size(), x1.size());
    // tensor label: I411
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

void Task454::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I411
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    sort_indices<0,1,1,1,-1,2>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task455::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c4 = b(1);
  const Index a2 = b(2);
  const Index c3 = b(3);
  // tensor label: I410
  std::unique_ptr<double[]> odata = out()->move_block(x0, c4, a2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c4, a2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c4, a2, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c4, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c4, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c4.size(), x1.size());
    // tensor label: I415
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

void Task456::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I415
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    sort_indices<0,1,1,1,1,4>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task457::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c1 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I391
  std::unique_ptr<double[]> odata = out()->move_block(x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c1, c3, a2), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, c4)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), c4.size());
    // tensor label: I418
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c4, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c4, a2, c1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c4.size(), a2.size(), c1.size());
    dgemm_("T", "N", c3.size(), x0.size()*a2.size()*c1.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c3.size());
  }
  sort_indices<1,3,0,2,1,1,1,1>(odata_sorted, odata, c3.size(), x0.size(), a2.size(), c1.size());
  out()->put_block(odata, x0, c1, c3, a2);
}

void Task458::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c4 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I418
  std::unique_ptr<double[]> odata = out()->move_block(x0, c4, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c4, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c4, a2, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c1, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a2, c1, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c1.size(), x1.size());
    // tensor label: I419
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

void Task459::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I419
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    sort_indices<0,1,1,1,1,4>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task460::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c4 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I418
  std::unique_ptr<double[]> odata = out()->move_block(x0, c4, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c4, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c4, a2, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c4, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c4, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c4.size(), x1.size());
    // tensor label: I427
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

void Task461::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I427
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    sort_indices<0,1,1,1,-1,2>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task462::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c1 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I391
  std::unique_ptr<double[]> odata = out()->move_block(x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c1, c3, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a2.size());
    // tensor label: I422
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a4, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a4, c1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a4.size(), c1.size());
    dgemm_("T", "N", a2.size(), x0.size()*c3.size()*c1.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<1,3,2,0,1,1,1,1>(odata_sorted, odata, a2.size(), x0.size(), c3.size(), c1.size());
  out()->put_block(odata, x0, c1, c3, a2);
}

void Task463::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a4 = b(2);
  const Index c1 = b(3);
  // tensor label: I422
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a4, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3, a4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3, a4, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, c1, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4, c1, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), c1.size(), x1.size());
    // tensor label: I423
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

void Task464::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I423
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    sort_indices<0,1,1,1,-1,4>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task465::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a4 = b(2);
  const Index c1 = b(3);
  // tensor label: I422
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a4, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3, a4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3, a4, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x1.size());
    // tensor label: I431
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

void Task466::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I431
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    sort_indices<0,1,1,1,1,2>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task467::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c1 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I391
  std::unique_ptr<double[]> odata = out()->move_block(x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c1, c3, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), x1.size());
    // tensor label: I434
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, a2, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, a2, c3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), a2.size(), c3.size());
    dgemm_("T", "N", c1.size(), x0.size()*a2.size()*c3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), x0.size(), a2.size(), c3.size());
  out()->put_block(odata, x0, c1, c3, a2);
}

void Task468::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c3 = b(3);
  // tensor label: I434
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, a2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a2, c3), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c3, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c3.size(), x2.size());
      // tensor label: I435
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x1, x0, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x1, x0, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x1.size(), x0.size(), x2.size());
      dgemm_("T", "N", a2.size()*c3.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a2.size()*c3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a2.size(), c3.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, a2, c3);
}

void Task469::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  // tensor label: I435
  std::unique_ptr<double[]> odata = out()->move_block(x3, x1, x0, x2);
  {
    // tensor label: Gamma22
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x0, x2);
    sort_indices<0,1,2,3,1,1,1,4>(i0data, odata, x3.size(), x1.size(), x0.size(), x2.size());
  }
  out()->put_block(odata, x3, x1, x0, x2);
}

void Task470::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c3 = b(3);
  // tensor label: I434
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, a2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a2, c3), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), x3.size(), x2.size());
      // tensor label: I443
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      dgemm_("T", "N", c3.size()*a2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c3.size()*a2.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), x0.size(), x1.size());
  out()->put_block(odata, x1, x0, a2, c3);
}

void Task471::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I443
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x0, x1);
  {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
    sort_indices<0,1,2,3,1,1,-1,4>(i0data, odata, x3.size(), x2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x3, x2, x0, x1);
}

void Task472::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c1 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I391
  std::unique_ptr<double[]> odata = out()->move_block(x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c1, c3, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), x1.size());
    // tensor label: I438
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, a2, c1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a2.size(), c1.size());
    dgemm_("T", "N", c3.size(), x0.size()*a2.size()*c1.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c3.size());
  }
  sort_indices<1,3,0,2,1,1,1,1>(odata_sorted, odata, c3.size(), x0.size(), a2.size(), c1.size());
  out()->put_block(odata, x0, c1, c3, a2);
}

void Task473::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I438
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a2, c1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c1.size(), x2.size());
      // tensor label: I439
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      dgemm_("T", "N", a2.size()*c1.size(), x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, a2, c1);
}

void Task474::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I439
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x0, x1);
  {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
    sort_indices<0,1,2,3,1,1,-1,4>(i0data, odata, x3.size(), x2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x3, x2, x0, x1);
}

void Task475::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I438
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a2, c1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x3.size(), x2.size());
      // tensor label: I447
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      dgemm_("T", "N", c1.size()*a2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, a2, c1);
}

void Task476::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I447
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x0, x1);
  {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x3.size(), x2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x3, x2, x0, x1);
}

void Task477::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c1 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I391
  std::unique_ptr<double[]> odata = out()->move_block(x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c1, c3, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, x1)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), x1.size());
      // tensor label: I450
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, c1, a4, c3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, c1, a4, c3, a2)]);
      sort_indices<3,1,0,2,4,5,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), c1.size(), a4.size(), c3.size(), a2.size());
      dgemm_("T", "N", 1, x0.size()*c1.size()*c3.size()*a2.size(), x1.size()*a4.size(),
             1.0, i0data_sorted, x1.size()*a4.size(), i1data_sorted, x1.size()*a4.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x0.size(), c1.size(), c3.size(), a2.size());
  out()->put_block(odata, x0, c1, c3, a2);
}

void Task478::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index a4 = b(3);
  const Index c3 = b(4);
  const Index a2 = b(5);
  // tensor label: I450
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, c1, a4, c3, a2), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
  // tensor label: I451
  std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
  dgemm_("T", "N", c1.size()*a4.size()*c3.size()*a2.size(), x0.size()*x1.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a4.size()*c3.size()*a2.size());
  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, c1, a4, c3, a2);
}

void Task479::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I451
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    sort_indices<0,1,1,1,-1,2>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task480::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index a4 = b(3);
  const Index c3 = b(4);
  const Index a2 = b(5);
  // tensor label: I450
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, c1, a4, c3, a2), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
  // tensor label: I455
  std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
  dgemm_("T", "N", c1.size()*a2.size()*c3.size()*a4.size(), x0.size()*x1.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a2.size()*c3.size()*a4.size());
  sort_indices<4,5,0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, c1, a4, c3, a2);
}

void Task481::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I455
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task482::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c1 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I391
  std::unique_ptr<double[]> odata = out()->move_block(x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c1, c3, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c1, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c1.size(), x1.size());
    // tensor label: I730
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
    dgemm_("T", "N", c3.size()*a2.size()*c1.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c3.size()*a2.size()*c1.size());
  }
  sort_indices<3,2,0,1,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), x0.size());
  out()->put_block(odata, x0, c1, c3, a2);
}

void Task483::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I730
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    dscal_(x0.size()*x1.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,1,1,1,4>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task484::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c1 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I391
  std::unique_ptr<double[]> odata = out()->move_block(x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c1, c3, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x1.size());
    // tensor label: I733
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
    dgemm_("T", "N", c1.size()*a2.size()*c3.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<3,0,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), x0.size());
  out()->put_block(odata, x0, c1, c3, a2);
}

void Task485::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I733
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    dscal_(x0.size()*x1.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,1,1,-1,2>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task486::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c1 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I391
  std::unique_ptr<double[]> odata = out()->move_block(x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c1, c3, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x1, c1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x1, c1, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), x1.size(), c1.size(), a2.size());
    // tensor label: I773
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
    dgemm_("T", "N", c3.size()*c1.size()*a2.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c3.size()*c1.size()*a2.size());
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, c3.size(), c1.size(), a2.size(), x0.size());
  out()->put_block(odata, x0, c1, c3, a2);
}

void Task487::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I773
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    sort_indices<0,1,1,1,2,1>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task488::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c1 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I391
  std::unique_ptr<double[]> odata = out()->move_block(x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c1, c3, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x1, c3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x1, c3, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x1.size(), c3.size(), a2.size());
    // tensor label: I776
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
    dgemm_("T", "N", c1.size()*c3.size()*a2.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size()*c3.size()*a2.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), a2.size(), x0.size());
  out()->put_block(odata, x0, c1, c3, a2);
}

void Task489::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I776
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task490::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a1 = b(1);
  const Index c2 = b(2);
  const Index x1 = b(3);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, x1)]);
  sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), x1.size());
  // tensor label: I457
  std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, c2, a1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, c2, a1)]);
  sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), c2.size(), a1.size());
  energy_ += ddot_(x1.size()*x0.size()*c2.size()*a1.size(), i0data_sorted, 1, i1data_sorted, 1);
}

void Task491::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I457
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), a1.size());
    // tensor label: I458
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x2, x0, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x2, x0, c2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), x2.size(), x0.size(), c2.size());
    dgemm_("T", "N", a1.size(), x1.size()*x0.size()*c2.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x1.size(), x0.size(), c2.size());
  out()->put_block(odata, x1, x0, c2, a1);
}

void Task492::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  // tensor label: I458
  std::unique_ptr<double[]> odata = out()->move_block(x1, x2, x0, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x2, x0, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x2, x0, c2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c2, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, c2, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), c2.size(), x3.size());
        // tensor label: I459
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x1, x3, x2, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x1, x3, x2, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x1.size(), x3.size(), x2.size(), x0.size());
        dgemm_("T", "N", c2.size(), x1.size()*x2.size()*x0.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, c2.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c2.size(), x1.size(), x2.size(), x0.size());
  out()->put_block(odata, x1, x2, x0, c2);
}

void Task493::Task_local::compute() {
  energy_ = 0.0;
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x1 = b(2);
  const Index x3 = b(3);
  const Index x2 = b(4);
  const Index x0 = b(5);
  // tensor label: I459
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, x1, x3, x2, x0);
  {
    // tensor label: Gamma28
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x1, x3, x2, x0);
    sort_indices<0,1,2,3,4,5,1,1,-1,4>(i0data, odata, x5.size(), x4.size(), x1.size(), x3.size(), x2.size(), x0.size());
  }
  out()->put_block(odata, x5, x4, x1, x3, x2, x0);
}

void Task494::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I457
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, c3)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), c3.size());
      // tensor label: I462
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x2, x0, c3, a1, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x2, x0, c3, a1, c2)]);
      sort_indices<1,3,0,2,4,5,0,1,1,1>(i1data, i1data_sorted, x1.size(), x2.size(), x0.size(), c3.size(), a1.size(), c2.size());
      dgemm_("T", "N", 1, x1.size()*x0.size()*a1.size()*c2.size(), x2.size()*c3.size(),
             1.0, i0data_sorted, x2.size()*c3.size(), i1data_sorted, x2.size()*c3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a1.size(), c2.size());
  out()->put_block(odata, x1, x0, c2, a1);
}

void Task495::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index c3 = b(3);
  const Index a1 = b(4);
  const Index c2 = b(5);
  // tensor label: I462
  std::unique_ptr<double[]> odata = out()->move_block(x1, x2, x0, c3, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x2, x0, c3, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x2, x0, c3, a1, c2), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a1, c2, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a1, c2, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), a1.size(), c2.size(), x3.size());
    // tensor label: I463
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x3, x2, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x3, x2, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), x3.size(), x2.size(), x0.size());
    dgemm_("T", "N", c3.size()*a1.size()*c2.size(), x1.size()*x2.size()*x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c3.size()*a1.size()*c2.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, c3.size(), a1.size(), c2.size(), x1.size(), x2.size(), x0.size());
  out()->put_block(odata, x1, x2, x0, c3, a1, c2);
}

void Task496::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x0 = b(3);
  // tensor label: I463
  std::unique_ptr<double[]> odata = out()->move_block(x1, x3, x2, x0);
  {
    // tensor label: Gamma29
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x3, x2, x0);
    sort_indices<0,1,2,3,1,1,1,4>(i0data, odata, x1.size(), x3.size(), x2.size(), x0.size());
  }
  out()->put_block(odata, x1, x3, x2, x0);
}

void Task497::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index c3 = b(3);
  const Index a1 = b(4);
  const Index c2 = b(5);
  // tensor label: I462
  std::unique_ptr<double[]> odata = out()->move_block(x1, x2, x0, c3, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x2, x0, c3, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x2, x0, c3, a1, c2), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c3, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c3, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c3.size(), x3.size());
    // tensor label: I467
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x3, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x3, x1, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
    dgemm_("T", "N", c2.size()*a1.size()*c3.size(), x2.size()*x1.size()*x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c2.size()*a1.size()*c3.size());
  }
  sort_indices<4,3,5,2,1,0,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), c3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x2, x0, c3, a1, c2);
}

void Task498::Task_local::compute() {
  energy_ = 0.0;
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I467
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, x1, x0);
  {
    // tensor label: Gamma7
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,4>(i0data, odata, x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x2, x3, x1, x0);
}

void Task499::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I457
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c2, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c2, x4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c2.size(), x4.size());
      // tensor label: I470
      std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x0, x1, x4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x0, x1, x4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x5.size(), x0.size(), x1.size(), x4.size());
      dgemm_("T", "N", a1.size()*c2.size(), x0.size()*x1.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), x0.size(), x1.size());
  out()->put_block(odata, x1, x0, c2, a1);
}

