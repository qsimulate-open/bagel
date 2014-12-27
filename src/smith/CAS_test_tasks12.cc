//
// BAGEL - Parallel electron correlation program.
// Filename: CAS_test_tasks12.cc
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


#include <src/smith/CAS_test_tasks12.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CAS_test;

void Task550::Task_local::compute() {
  energy_ = 0.0;
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I524
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, x1, x0);
  {
    // tensor label: Gamma34
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,4>(i0data, odata, x5.size(), x4.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x5, x4, x1, x0);
}

void Task551::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I511
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, c3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), c3.size());
    // tensor label: I527
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

void Task552::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c3 = b(3);
  // tensor label: I527
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, a2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a2, c3), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c3, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c3.size(), x2.size());
      // tensor label: I528
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

void Task553::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I528
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,1,4>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task554::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c3 = b(3);
  // tensor label: I527
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, a2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a2, c3), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), x3.size(), x2.size());
      // tensor label: I539
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

void Task555::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I539
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task556::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I511
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), a2.size());
    // tensor label: I531
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

void Task557::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a3 = b(2);
  const Index c1 = b(3);
  // tensor label: I531
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a3, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, a3, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a3, c1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c1.size(), x2.size());
      // tensor label: I532
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

void Task558::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I532
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,4>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task559::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a3 = b(2);
  const Index c1 = b(3);
  // tensor label: I531
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a3, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, a3, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a3, c1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), x3.size(), x2.size());
      // tensor label: I543
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

void Task560::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I543
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task561::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I511
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x5, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x5, x4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x5.size(), x4.size());
      // tensor label: I535
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

void Task562::Task_local::compute() {
  energy_ = 0.0;
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I535
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, x1, x0);
  {
    // tensor label: Gamma34
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x1, x0);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x5.size(), x4.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x5, x4, x1, x0);
}

void Task563::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I511
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), x2.size());
    // tensor label: I546
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

void Task564::Task_local::compute() {
  energy_ = 0.0;
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I546
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
        // tensor label: I547
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

void Task565::Task_local::compute() {
  energy_ = 0.0;
  const Index x5 = b(0);
  const Index x2 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I547
  std::unique_ptr<double[]> odata = out()->move_block(x5, x2, x4, x3, x1, x0);
  {
    // tensor label: Gamma51
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2, x4, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,1,1,-1,4>(i0data, odata, x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x5, x2, x4, x3, x1, x0);
}

void Task566::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I511
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, c3)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), c3.size());
      // tensor label: I550
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

void Task567::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a4 = b(3);
  const Index c3 = b(4);
  const Index a2 = b(5);
  // tensor label: I550
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a4, c3, a2), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
  // tensor label: I551
  std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
  dgemm_("T", "N", c1.size()*a4.size()*c3.size()*a2.size(), x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a4.size()*c3.size()*a2.size());
  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c1, a4, c3, a2);
}

void Task568::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I551
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task569::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a4 = b(3);
  const Index c3 = b(4);
  const Index a2 = b(5);
  // tensor label: I550
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a4, c3, a2), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
  // tensor label: I555
  std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
  dgemm_("T", "N", c1.size()*a2.size()*c3.size()*a4.size(), x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a2.size()*c3.size()*a4.size());
  sort_indices<4,5,0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c1, a4, c3, a2);
}

void Task570::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I555
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,2,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task571::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I511
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, x2)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), x2.size());
      // tensor label: I558
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

void Task572::Task_local::compute() {
  energy_ = 0.0;
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a3 = b(3);
  const Index c1 = b(4);
  const Index a2 = b(5);
  // tensor label: I558
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, a3, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, a3, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, a3, c1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c1, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c1.size(), a2.size());
    // tensor label: I559
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

void Task573::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I559
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task574::Task_local::compute() {
  energy_ = 0.0;
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a3 = b(3);
  const Index c1 = b(4);
  const Index a2 = b(5);
  // tensor label: I558
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, a3, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, a3, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, a3, c1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c1, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c1.size(), a3.size());
    // tensor label: I563
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
    dgemm_("T", "N", a2.size()*c1.size()*a3.size(), x2.size()*x1.size()*x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a2.size()*c1.size()*a3.size());
  }
  sort_indices<3,4,5,2,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), a3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x1, x0, a3, c1, a2);
}

void Task575::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I563
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,4>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task576::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I511
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c1.size(), x2.size());
      // tensor label: I742
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", a2.size()*c1.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c1, a2);
}

void Task577::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I742
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    dscal_(x3.size()*x2.size()*x1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,1,1,1,4>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task578::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I511
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x3.size(), x2.size());
      // tensor label: I745
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c1, a2);
}

void Task579::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I745
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    dscal_(x3.size()*x2.size()*x1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task580::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I511
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x3.size(), x2.size());
      // tensor label: I791
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c1, a2);
}

void Task581::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I791
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task582::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I511
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x3, x2, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x3, x2, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x3.size(), x2.size(), a2.size());
      // tensor label: I794
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x3, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x3, x1, x0)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), x1.size()*x0.size(), x2.size()*x3.size(),
             1.0, i0data_sorted, x2.size()*x3.size(), i1data_sorted, x2.size()*x3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c1, a2);
}

void Task583::Task_local::compute() {
  energy_ = 0.0;
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I794
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, x1, x0);
  {
    // tensor label: Gamma7
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x2, x3, x1, x0);
}

void Task584::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I511
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c1.size(), x2.size());
      // tensor label: I797
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", a2.size()*c1.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c1, a2);
}

void Task585::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I797
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task586::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I511
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, c1, a2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), c1.size(), a2.size());
      // tensor label: I800
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c1, a2);
}

void Task587::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I800
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task588::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I511
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  // tensor label: h1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size());
  // tensor label: I828
  std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
  dgemm_("T", "N", c1.size()*a2.size(), x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a2.size());
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c1, a2);
}

void Task589::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I828
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,2,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task590::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a1 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
  sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
  // tensor label: I565
  std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x2, x1, a1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x2, x1, a1)]);
  sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x2.size(), x1.size(), a1.size());
  energy_ += ddot_(x0.size()*x2.size()*x1.size()*a1.size(), i0data_sorted, 1, i1data_sorted, 1);
}

void Task591::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  // tensor label: I565
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, c2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, c2)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), c2.size());
      // tensor label: I566
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x3, x2, x1, a1, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x3, x2, x1, a1, c2)]);
      sort_indices<1,5,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, x0.size(), x3.size(), x2.size(), x1.size(), a1.size(), c2.size());
      dgemm_("T", "N", 1, x0.size()*x2.size()*x1.size()*a1.size(), x3.size()*c2.size(),
             1.0, i0data_sorted, x3.size()*c2.size(), i1data_sorted, x3.size()*c2.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->put_block(odata, x0, x2, x1, a1);
}

void Task592::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a1 = b(4);
  const Index c2 = b(5);
  // tensor label: I566
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3, x2, x1, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3, x2, x1, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3, x2, x1, a1, c2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c2, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c2, x4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c2.size(), x4.size());
      // tensor label: I567
      std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x0, x3, x4, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x0, x3, x4, x2, x1)]);
      sort_indices<0,3,1,2,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x0.size(), x3.size(), x4.size(), x2.size(), x1.size());
      dgemm_("T", "N", a1.size()*c2.size(), x0.size()*x3.size()*x2.size()*x1.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<2,3,4,5,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), x0.size(), x3.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x3, x2, x1, a1, c2);
}

void Task593::Task_local::compute() {
  energy_ = 0.0;
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);
  // tensor label: I567
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x3, x4, x2, x1);
  {
    // tensor label: Gamma56
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x3, x4, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,1,4>(i0data, odata, x5.size(), x0.size(), x3.size(), x4.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x5, x0, x3, x4, x2, x1);
}

void Task594::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a1 = b(4);
  const Index c2 = b(5);
  // tensor label: I566
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3, x2, x1, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3, x2, x1, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3, x2, x1, a1, c2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x5, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x5, x4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x5.size(), x4.size());
      // tensor label: I571
      std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x3, x0, x2, x1)]);
      sort_indices<0,1,2,3,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
      dgemm_("T", "N", c2.size()*a1.size(), x3.size()*x0.size()*x2.size()*x1.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<3,2,4,5,1,0,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x3.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x3, x2, x1, a1, c2);
}

void Task595::Task_local::compute() {
  energy_ = 0.0;
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);
  // tensor label: I571
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, x3, x0, x2, x1);
  {
    // tensor label: Gamma57
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,-1,4>(i0data, odata, x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x5, x4, x3, x0, x2, x1);
}

void Task596::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  // tensor label: I565
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a1, x6, x5);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a1, x6, x5)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x7.size(), a1.size(), x6.size(), x5.size());
        // tensor label: I574
        std::unique_ptr<double[]> i1data = in(1)->get_block(x7, x0, x6, x5, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x7, x0, x6, x5, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i1data, i1data_sorted, x7.size(), x0.size(), x6.size(), x5.size(), x2.size(), x1.size());
        dgemm_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a1);
}

void Task597::Task_local::compute() {
  energy_ = 0.0;
  const Index x7 = b(0);
  const Index x0 = b(1);
  const Index x6 = b(2);
  const Index x5 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);
  // tensor label: I574
  std::unique_ptr<double[]> odata = out()->move_block(x7, x0, x6, x5, x2, x1);
  {
    // tensor label: Gamma58
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x5, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,1,4>(i0data, odata, x7.size(), x0.size(), x6.size(), x5.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x7, x0, x6, x5, x2, x1);
}

void Task598::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  // tensor label: I565
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a2, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a2.size(), a1.size());
    // tensor label: I577
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x2, x1, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x2, x1, a2)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), x2.size(), x1.size(), a2.size());
    dgemm_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a1);
}

void Task599::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  // tensor label: I577
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, x4, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), x4.size(), x3.size());
        // tensor label: I578
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        dgemm_("T", "N", a2.size(), x0.size()*x2.size()*x1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a2);
}

