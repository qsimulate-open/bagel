//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks11.cc
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


#include <src/smith/CASPT2_tasks11.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task500::Task_local::compute() {
  target_ = 0.0;
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I652
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task501::Task_local::compute() {
  target_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I487
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, x1)]);
      sort_indices<2,1,0,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), x1.size());
      // tensor label: I504
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x3, c2, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x3, c2, a1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, x2.size(), x3.size(), c2.size(), a1.size());
      dgemm_("T", "N", x1.size()*x0.size(), x2.size()*x3.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size());
  out()->put_block(odata, x1, x0, x3, x2);
}

void Task502::Task_local::compute() {
  target_ = 0.0;
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I504
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, c2, a1), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a1.size());
    // tensor label: I505
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x2)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), x2.size());
    dgemm_("T", "N", x3.size()*c2.size()*a1.size(), x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x3.size()*c2.size()*a1.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x3.size(), c2.size(), a1.size(), x2.size());
  out()->put_block(odata, x2, x3, c2, a1);
}

void Task503::Task_local::compute() {
  target_ = 0.0;
  const Index a3 = b(0);
  const Index x2 = b(1);
  // tensor label: I505
  std::unique_ptr<double[]> odata = out()->move_block(a3, x2);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x2);
    sort_indices<0,1,1,1,-1,4>(i0data, odata, a3.size(), x2.size());
  }
  out()->put_block(odata, a3, x2);
}

void Task504::Task_local::compute() {
  target_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I487
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c1, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c1.size(), x2.size());
      // tensor label: I531
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, c1, a3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, c1, a3)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), c1.size(), a3.size());
      dgemm_("T", "N", x3.size()*x2.size(), x1.size()*x0.size(), c1.size()*a3.size(),
             1.0, i0data_sorted, c1.size()*a3.size(), i1data_sorted, c1.size()*a3.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, x3, x2);
}

void Task505::Task_local::compute() {
  target_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a3 = b(3);
  // tensor label: I531
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a3), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, a2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a3.size(), a2.size());
    // tensor label: I532
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, a2, c1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), a2.size(), c1.size());
    dgemm_("T", "N", a3.size(), x1.size()*x0.size()*c1.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, a3.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a3.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, x1, x0, c1, a3);
}

void Task506::Task_local::compute() {
  target_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I532
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
    sort_indices<3,2,1,0,1,1,-1,4>(i0data, odata, c1.size(), a2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x1, x0, a2, c1);
}

void Task507::Task_local::compute() {
  target_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I487
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, x3, x2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), x3.size(), x2.size());
      // tensor label: I542
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, c1, a3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, c1, a3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), c1.size(), a3.size());
      dgemm_("T", "N", x3.size()*x2.size(), x1.size()*x0.size(), c1.size()*a3.size(),
             1.0, i0data_sorted, c1.size()*a3.size(), i1data_sorted, c1.size()*a3.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, x3, x2);
}

void Task508::Task_local::compute() {
  target_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a3 = b(3);
  // tensor label: I542
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a3), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, a2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a3.size(), a2.size());
    // tensor label: I543
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, a2, c1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), a2.size(), c1.size());
    dgemm_("T", "N", a3.size(), x1.size()*x0.size()*c1.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, a3.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a3.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, x1, x0, c1, a3);
}

void Task509::Task_local::compute() {
  target_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I543
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, c1.size(), a2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x1, x0, a2, c1);
}

void Task510::Task_local::compute() {
  target_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I487
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x0, x1)]);
      sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x0.size(), x1.size());
      // tensor label: I558
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x3, c1, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x3, c1, a2)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, x2.size(), x3.size(), c1.size(), a2.size());
      dgemm_("T", "N", x1.size()*x0.size(), x2.size()*x3.size(), c1.size()*a2.size(),
             1.0, i0data_sorted, c1.size()*a2.size(), i1data_sorted, c1.size()*a2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size());
  out()->put_block(odata, x1, x0, x3, x2);
}

void Task511::Task_local::compute() {
  target_ = 0.0;
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I558
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x3, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, c1, a2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c1, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c1.size(), a2.size());
    // tensor label: I559
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x2)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), x2.size());
    dgemm_("T", "N", x3.size()*c1.size()*a2.size(), x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x3.size()*c1.size()*a2.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x3.size(), c1.size(), a2.size(), x2.size());
  out()->put_block(odata, x2, x3, c1, a2);
}

void Task512::Task_local::compute() {
  target_ = 0.0;
  const Index a3 = b(0);
  const Index x2 = b(1);
  // tensor label: I559
  std::unique_ptr<double[]> odata = out()->move_block(a3, x2);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x2);
    sort_indices<0,1,1,1,1,2>(i0data, odata, a3.size(), x2.size());
  }
  out()->put_block(odata, a3, x2);
}

void Task513::Task_local::compute() {
  target_ = 0.0;
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I558
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x3, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, c1, a2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c1, a3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c1.size(), a3.size());
    // tensor label: I563
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x2)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), x2.size());
    dgemm_("T", "N", x3.size()*a2.size()*c1.size(), x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x3.size()*a2.size()*c1.size());
  }
  sort_indices<3,0,2,1,1,1,1,1>(odata_sorted, odata, x3.size(), a2.size(), c1.size(), x2.size());
  out()->put_block(odata, x2, x3, c1, a2);
}

void Task514::Task_local::compute() {
  target_ = 0.0;
  const Index a3 = b(0);
  const Index x2 = b(1);
  // tensor label: I563
  std::unique_ptr<double[]> odata = out()->move_block(a3, x2);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x2);
    sort_indices<0,1,1,1,-1,4>(i0data, odata, a3.size(), x2.size());
  }
  out()->put_block(odata, a3, x2);
}

void Task515::Task_local::compute() {
  target_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I487
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), x2.size());
      // tensor label: I643
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), x0.size(), x1.size());
      dgemm_("T", "N", x3.size()*x2.size(), x0.size()*x1.size(), a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<3,2,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x0.size(), x1.size());
  out()->put_block(odata, x1, x0, x3, x2);
}

void Task516::Task_local::compute() {
  target_ = 0.0;
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I643
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, x1), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size());
    // tensor label: I644
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x0)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x0.size());
    dgemm_("T", "N", x1.size(), a3.size()*c2.size()*x0.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), a3.size(), c2.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, x1);
}

void Task517::Task_local::compute() {
  target_ = 0.0;
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I644
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,-1,4>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task518::Task_local::compute() {
  target_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I487
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x3, x2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x3.size(), x2.size());
      // tensor label: I655
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x0.size(), x1.size());
      dgemm_("T", "N", x3.size()*x2.size(), x0.size()*x1.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<3,2,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x0.size(), x1.size());
  out()->put_block(odata, x1, x0, x3, x2);
}

void Task519::Task_local::compute() {
  target_ = 0.0;
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I655
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, x0, x1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
    dscal_(x1.size()*c2.size()*a1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<2,1,0,3,1,1,1,4>(i0data, odata, x0.size(), a1.size(), c2.size(), x1.size());
  }
  out()->put_block(odata, c2, a1, x0, x1);
}

void Task520::Task_local::compute() {
  target_ = 0.0;
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I655
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x0, x1), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size());
    // tensor label: I656
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x0)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x0.size());
    dgemm_("T", "N", x1.size(), c2.size()*a1.size()*x0.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), c2.size(), a1.size(), x0.size());
  out()->put_block(odata, c2, a1, x0, x1);
}

void Task521::Task_local::compute() {
  target_ = 0.0;
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I656
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,-1,4>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task522::Task_local::compute() {
  target_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I487
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c1, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c1.size(), x2.size());
      // tensor label: I742
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

void Task523::Task_local::compute() {
  target_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I742
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
    dscal_(x1.size()*x0.size()*a2.size()*c1.size(), e0_, i0data.get(), 1);
    sort_indices<3,2,1,0,1,1,1,4>(i0data, odata, c1.size(), a2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x1, x0, a2, c1);
}

void Task524::Task_local::compute() {
  target_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I487
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x3, x2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x3.size(), x2.size());
      // tensor label: I745
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

void Task525::Task_local::compute() {
  target_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I745
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
    dscal_(x1.size()*x0.size()*a2.size()*c1.size(), e0_, i0data.get(), 1);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, c1.size(), a2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x1, x0, a2, c1);
}

void Task526::Task_local::compute() {
  target_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I487
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x3, x2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x3.size(), x2.size());
      // tensor label: I779
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

void Task527::Task_local::compute() {
  target_ = 0.0;
  const Index x1 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I779
  std::unique_ptr<double[]> odata = out()->move_block(x1, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, x0.size(), a1.size(), c2.size(), x1.size());
  }
  out()->put_block(odata, x1, c2, a1, x0);
}

void Task528::Task_local::compute() {
  target_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I487
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, c2, a1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), c2.size(), a1.size());
      // tensor label: I788
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

void Task529::Task_local::compute() {
  target_ = 0.0;
  const Index x1 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I788
  std::unique_ptr<double[]> odata = out()->move_block(x1, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, x0.size(), a1.size(), c2.size(), x1.size());
  }
  out()->put_block(odata, x1, c2, a1, x0);
}

void Task530::Task_local::compute() {
  target_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I487
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x3, x2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x3.size(), x2.size());
      // tensor label: I791
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

void Task531::Task_local::compute() {
  target_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I791
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, c1.size(), a2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x1, x0, a2, c1);
}

void Task532::Task_local::compute() {
  target_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I487
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c1, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c1.size(), x2.size());
      // tensor label: I797
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

void Task533::Task_local::compute() {
  target_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I797
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, c1.size(), a2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x1, x0, a2, c1);
}

void Task534::Task_local::compute() {
  target_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I487
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, c1, a2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), c1.size(), a2.size());
      // tensor label: I800
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

void Task535::Task_local::compute() {
  target_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I800
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, c1.size(), a2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x1, x0, a2, c1);
}

void Task536::Task_local::compute() {
  target_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: Gamma38
  std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
  // tensor label: I495
  std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
  target_ += ddot_(x1.size()*x0.size(), i0data_sorted, 1, i1data_sorted, 1);
}

void Task537::Task_local::compute() {
  target_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I495
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, x1)]);
      sort_indices<2,1,0,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), x1.size());
      // tensor label: I496
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

void Task538::Task_local::compute() {
  target_ = 0.0;
  const Index c2 = b(0);
  const Index a1 = b(1);
  // tensor label: I496
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a4, c3, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a4, c3, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a4.size(), c3.size(), a1.size());
      // tensor label: I497
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

void Task539::Task_local::compute() {
  target_ = 0.0;
  const Index a4 = b(0);
  const Index c3 = b(1);
  // tensor label: I497
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3);
    sort_indices<0,1,1,1,1,2>(i0data, odata, a4.size(), c3.size());
  }
  out()->put_block(odata, a4, c3);
}

void Task540::Task_local::compute() {
  target_ = 0.0;
  const Index c2 = b(0);
  const Index a1 = b(1);
  // tensor label: I496
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c3, a4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c3.size(), a4.size());
      // tensor label: I501
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

void Task541::Task_local::compute() {
  target_ = 0.0;
  const Index a4 = b(0);
  const Index c3 = b(1);
  // tensor label: I501
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, a4.size(), c3.size());
  }
  out()->put_block(odata, a4, c3);
}

void Task542::Task_local::compute() {
  target_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I495
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x0, x1)]);
      sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x0.size(), x1.size());
      // tensor label: I550
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

void Task543::Task_local::compute() {
  target_ = 0.0;
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: I550
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
      // tensor label: I551
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

void Task544::Task_local::compute() {
  target_ = 0.0;
  const Index a4 = b(0);
  const Index c3 = b(1);
  // tensor label: I551
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, a4.size(), c3.size());
  }
  out()->put_block(odata, a4, c3);
}

void Task545::Task_local::compute() {
  target_ = 0.0;
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: I550
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
      // tensor label: I555
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

void Task546::Task_local::compute() {
  target_ = 0.0;
  const Index a4 = b(0);
  const Index c3 = b(1);
  // tensor label: I555
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3);
    sort_indices<0,1,1,1,2,1>(i0data, odata, a4.size(), c3.size());
  }
  out()->put_block(odata, a4, c3);
}

void Task547::Task_local::compute() {
  target_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I495
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c1, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c1.size(), x0.size());
      // tensor label: I601
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

void Task548::Task_local::compute() {
  target_ = 0.0;
  const Index a4 = b(0);
  const Index c1 = b(1);
  // tensor label: I601
  std::unique_ptr<double[]> odata = out()->move_block(a4, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
      sort_indices<2,1,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
      // tensor label: I602
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size());
      dgemm_("T", "N", a4.size()*c1.size(), 1, c3.size()*a2.size(),
             1.0, i0data_sorted, c3.size()*a2.size(), i1data_sorted, c3.size()*a2.size(),
             1.0, odata_sorted, a4.size()*c1.size());
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size());
  out()->put_block(odata, a4, c1);
}

void Task549::Task_local::compute() {
  target_ = 0.0;
  const Index c3 = b(0);
  const Index a2 = b(1);
  // tensor label: I602
  std::unique_ptr<double[]> odata = out()->move_block(c3, a2);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2);
    sort_indices<0,1,1,1,1,2>(i0data, odata, c3.size(), a2.size());
  }
  out()->put_block(odata, c3, a2);
}

