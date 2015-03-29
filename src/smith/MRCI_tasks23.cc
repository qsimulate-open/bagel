//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks23.cc
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

#include <src/smith/MRCI_tasks23.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task1100::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I219
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& a5 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a5, c4, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a5, c4, a3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a5.size(), c4.size(), a3.size());
      // tensor label: I1591
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, a5, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, a5, c4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), a5.size(), c4.size());
      dgemm_("T", "N", x1.size()*a3.size(), c2.size()*a1.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, x1.size()*a3.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x1.size(), a3.size(), c2.size(), a1.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task1101::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index a5 = b(2);
  const Index c4 = b(3);
  // tensor label: I1591
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, a5, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, a5, c4);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c2.size(), a1.size(), a5.size(), c4.size());
  }
  out()->put_block(odata, c2, a1, a5, c4);
}

void Task1102::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I219
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c4, a5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c4, a5)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c4.size(), a5.size());
      // tensor label: I1594
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, a5, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, a5, c4)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), a5.size(), c4.size());
      dgemm_("T", "N", x1.size()*a3.size(), c2.size()*a1.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, x1.size()*a3.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x1.size(), a3.size(), c2.size(), a1.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task1103::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index a5 = b(2);
  const Index c4 = b(3);
  // tensor label: I1594
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, a5, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, a5, c4);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c2.size(), a1.size(), a5.size(), c4.size());
  }
  out()->put_block(odata, c2, a1, a5, c4);
}

void Task1104::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I219
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& a5 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a5, c4, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a5, c4, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a5.size(), c4.size(), a1.size());
      // tensor label: I1597
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a3, a5, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a3, a5, c4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), a3.size(), a5.size(), c4.size());
      dgemm_("T", "N", x1.size()*a1.size(), c2.size()*a3.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, x1.size()*a1.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, x1.size(), a1.size(), c2.size(), a3.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task1105::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index a5 = b(2);
  const Index c4 = b(3);
  // tensor label: I1597
  std::unique_ptr<double[]> odata = out()->move_block(c2, a3, a5, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, a5, c4);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c2.size(), a3.size(), a5.size(), c4.size());
  }
  out()->put_block(odata, c2, a3, a5, c4);
}

void Task1106::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I219
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c4, a5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c4, a5)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c4.size(), a5.size());
      // tensor label: I1600
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a3, a5, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a3, a5, c4)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), a3.size(), a5.size(), c4.size());
      dgemm_("T", "N", x1.size()*a1.size(), c2.size()*a3.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, x1.size()*a1.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, x1.size(), a1.size(), c2.size(), a3.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task1107::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index a5 = b(2);
  const Index c4 = b(3);
  // tensor label: I1600
  std::unique_ptr<double[]> odata = out()->move_block(c2, a3, a5, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, a5, c4);
    sort_indices<0,1,2,3,1,1,4,1>(i0data, odata, c2.size(), a3.size(), a5.size(), c4.size());
  }
  out()->put_block(odata, c2, a3, a5, c4);
}

void Task1108::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I219
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& a5 : *range_[2]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a5, c2, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a5, c2, a4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a5.size(), c2.size(), a4.size());
      // tensor label: I1603
      std::unique_ptr<double[]> i1data = in(1)->get_block(a5, a1, a4, a3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, a1, a4, a3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, a5.size(), a1.size(), a4.size(), a3.size());
      dgemm_("T", "N", x1.size()*c2.size(), a1.size()*a3.size(), a5.size()*a4.size(),
             1.0, i0data_sorted, a5.size()*a4.size(), i1data_sorted, a5.size()*a4.size(),
             1.0, odata_sorted, x1.size()*c2.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x1.size(), c2.size(), a1.size(), a3.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task1109::Task_local::compute() {
  const Index a5 = b(0);
  const Index a1 = b(1);
  const Index a4 = b(2);
  const Index a3 = b(3);
  // tensor label: I1603
  std::unique_ptr<double[]> odata = out()->move_block(a5, a1, a4, a3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, a1, a4, a3);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, a5.size(), a1.size(), a4.size(), a3.size());
  }
  out()->put_block(odata, a5, a1, a4, a3);
}

void Task1110::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I219
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a5)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a5.size());
      // tensor label: I1606
      std::unique_ptr<double[]> i1data = in(1)->get_block(a5, a1, a4, a3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, a1, a4, a3)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a5.size(), a1.size(), a4.size(), a3.size());
      dgemm_("T", "N", x1.size()*c2.size(), a1.size()*a3.size(), a5.size()*a4.size(),
             1.0, i0data_sorted, a5.size()*a4.size(), i1data_sorted, a5.size()*a4.size(),
             1.0, odata_sorted, x1.size()*c2.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x1.size(), c2.size(), a1.size(), a3.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task1111::Task_local::compute() {
  const Index a5 = b(0);
  const Index a1 = b(1);
  const Index a4 = b(2);
  const Index a3 = b(3);
  // tensor label: I1606
  std::unique_ptr<double[]> odata = out()->move_block(a5, a1, a4, a3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, a1, a4, a3);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, a5.size(), a1.size(), a4.size(), a3.size());
  }
  out()->put_block(odata, a5, a1, a4, a3);
}

void Task1112::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), x1.size());
    // tensor label: I237
    std::unique_ptr<double[]> i1data = in(1)->get_block(a1, a3, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, a3, x0, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a1.size(), a3.size(), x0.size(), x1.size());
    dgemm_("T", "N", c2.size(), a1.size()*a3.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c2.size());
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1113::Task_local::compute() {
  const Index a1 = b(0);
  const Index a3 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I237
  std::unique_ptr<double[]> odata = out()->move_block(a1, a3, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, a3, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, a3, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma51
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I238
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, x2, a3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, x2, a3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), x2.size(), a3.size());
      dgemm_("T", "N", x0.size()*x1.size(), a1.size()*a3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), a3.size());
  out()->put_block(odata, a1, a3, x0, x1);
}

void Task1114::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index x2 = b(2);
  const Index a3 = b(3);
  // tensor label: I238
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1, x2, a3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x3.size(), a1.size(), x2.size(), a3.size());
  }
  out()->put_block(odata, x3, a1, x2, a3);
}

void Task1115::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a1, x1, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a1, x1, a3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x2.size(), a1.size(), x1.size(), a3.size());
      // tensor label: I1359
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x1, x2, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x1, x2, x0)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), x1.size(), x2.size(), x0.size());
      dgemm_("T", "N", a1.size()*a3.size(), c2.size()*x0.size(), x1.size()*x2.size(),
             1.0, i0data_sorted, x1.size()*x2.size(), i1data_sorted, x1.size()*x2.size(),
             1.0, odata_sorted, a1.size()*a3.size());
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), c2.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1116::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index x0 = b(3);
  // tensor label: I1359
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, x2, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x2, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x2, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma24
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x1, x3, x2, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x1, x3, x2, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x1.size(), x3.size(), x2.size(), x0.size());
        // tensor label: I1360
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, c2, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, c2, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), c2.size(), x3.size());
        dgemm_("T", "N", x1.size()*x2.size()*x0.size(), c2.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x1.size()*x2.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x2.size(), x0.size(), c2.size());
  out()->put_block(odata, c2, x1, x2, x0);
}

void Task1117::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index c2 = b(2);
  const Index x3 = b(3);
  // tensor label: I1360
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, c2, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c2, x3);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x5.size(), x4.size(), c2.size(), x3.size());
  }
  out()->put_block(odata, x5, x4, c2, x3);
}

void Task1118::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a3, c2, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a3, c2, x3)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, c4.size(), a3.size(), c2.size(), x3.size());
      // tensor label: I1362
      std::unique_ptr<double[]> i1data = in(1)->get_block(a1, c4, x3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, c4, x3, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, a1.size(), c4.size(), x3.size(), x0.size());
      dgemm_("T", "N", a3.size()*c2.size(), a1.size()*x0.size(), c4.size()*x3.size(),
             1.0, i0data_sorted, c4.size()*x3.size(), i1data_sorted, c4.size()*x3.size(),
             1.0, odata_sorted, a3.size()*c2.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1119::Task_local::compute() {
  const Index a1 = b(0);
  const Index c4 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1362
  std::unique_ptr<double[]> odata = out()->move_block(a1, c4, x3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, c4, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, c4, x3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma25
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x3, x2, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x3, x2, x0)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), x3.size(), x2.size(), x0.size());
      // tensor label: I1363
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a1, x1, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a1, x1, c4)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), a1.size(), x1.size(), c4.size());
      dgemm_("T", "N", x3.size()*x0.size(), a1.size()*c4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a1.size(), c4.size());
  out()->put_block(odata, a1, c4, x3, x0);
}

void Task1120::Task_local::compute() {
  const Index x2 = b(0);
  const Index a1 = b(1);
  const Index x1 = b(2);
  const Index c4 = b(3);
  // tensor label: I1363
  std::unique_ptr<double[]> odata = out()->move_block(x2, a1, x1, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a1, x1, c4);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x2.size(), a1.size(), x1.size(), c4.size());
  }
  out()->put_block(odata, x2, a1, x1, c4);
}

void Task1121::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a1, c2, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a1, c2, x3)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, c4.size(), a1.size(), c2.size(), x3.size());
      // tensor label: I1365
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c4, x3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c4, x3, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c4.size(), x3.size(), x0.size());
      dgemm_("T", "N", a1.size()*c2.size(), a3.size()*x0.size(), c4.size()*x3.size(),
             1.0, i0data_sorted, c4.size()*x3.size(), i1data_sorted, c4.size()*x3.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<2,1,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1122::Task_local::compute() {
  const Index a3 = b(0);
  const Index c4 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1365
  std::unique_ptr<double[]> odata = out()->move_block(a3, c4, x3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c4, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c4, x3, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma5
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x3, x1, x0)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
      // tensor label: I1366
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a3, x1, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a3, x1, c4)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), a3.size(), x1.size(), c4.size());
      dgemm_("T", "N", x3.size()*x0.size(), a3.size()*c4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a3.size(), c4.size());
  out()->put_block(odata, a3, c4, x3, x0);
}

void Task1123::Task_local::compute() {
  const Index x2 = b(0);
  const Index a3 = b(1);
  const Index x1 = b(2);
  const Index c4 = b(3);
  // tensor label: I1366
  std::unique_ptr<double[]> odata = out()->move_block(x2, a3, x1, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a3, x1, c4);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x2.size(), a3.size(), x1.size(), c4.size());
  }
  out()->put_block(odata, x2, a3, x1, c4);
}

void Task1124::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c4, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c4, x3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), x3.size());
      // tensor label: I1368
      std::unique_ptr<double[]> i1data = in(1)->get_block(a1, c4, x3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, c4, x3, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, a1.size(), c4.size(), x3.size(), x0.size());
      dgemm_("T", "N", c2.size()*a3.size(), a1.size()*x0.size(), c4.size()*x3.size(),
             1.0, i0data_sorted, c4.size()*x3.size(), i1data_sorted, c4.size()*x3.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), a1.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1125::Task_local::compute() {
  const Index a1 = b(0);
  const Index c4 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1368
  std::unique_ptr<double[]> odata = out()->move_block(a1, c4, x3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, c4, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, c4, x3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma25
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x3, x2, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x3, x2, x0)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), x3.size(), x2.size(), x0.size());
      // tensor label: I1369
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a1, x1, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a1, x1, c4)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), a1.size(), x1.size(), c4.size());
      dgemm_("T", "N", x3.size()*x0.size(), a1.size()*c4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a1.size(), c4.size());
  out()->put_block(odata, a1, c4, x3, x0);
}

void Task1126::Task_local::compute() {
  const Index x2 = b(0);
  const Index a1 = b(1);
  const Index x1 = b(2);
  const Index c4 = b(3);
  // tensor label: I1369
  std::unique_ptr<double[]> odata = out()->move_block(x2, a1, x1, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a1, x1, c4);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x2.size(), a1.size(), x1.size(), c4.size());
  }
  out()->put_block(odata, x2, a1, x1, c4);
}

void Task1127::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c4, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c4, x3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c4.size(), x3.size());
      // tensor label: I1371
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c4, x3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c4, x3, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c4.size(), x3.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), a3.size()*x0.size(), c4.size()*x3.size(),
             1.0, i0data_sorted, c4.size()*x3.size(), i1data_sorted, c4.size()*x3.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1128::Task_local::compute() {
  const Index a3 = b(0);
  const Index c4 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1371
  std::unique_ptr<double[]> odata = out()->move_block(a3, c4, x3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c4, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c4, x3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma25
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x3, x2, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x3, x2, x0)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), x3.size(), x2.size(), x0.size());
      // tensor label: I1372
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a3, x1, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a3, x1, c4)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), a3.size(), x1.size(), c4.size());
      dgemm_("T", "N", x3.size()*x0.size(), a3.size()*c4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a3.size(), c4.size());
  out()->put_block(odata, a3, c4, x3, x0);
}

void Task1129::Task_local::compute() {
  const Index x2 = b(0);
  const Index a3 = b(1);
  const Index x1 = b(2);
  const Index c4 = b(3);
  // tensor label: I1372
  std::unique_ptr<double[]> odata = out()->move_block(x2, a3, x1, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a3, x1, c4);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x2.size(), a3.size(), x1.size(), c4.size());
  }
  out()->put_block(odata, x2, a3, x1, c4);
}

void Task1130::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, c2, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a3, c2, x4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a3.size(), c2.size(), x4.size());
      // tensor label: I1374
      std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x5, x4, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x5, x4, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, a1.size(), x5.size(), x4.size(), x0.size());
      dgemm_("T", "N", a3.size()*c2.size(), a1.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, a3.size()*c2.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1131::Task_local::compute() {
  const Index a1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  // tensor label: I1374
  std::unique_ptr<double[]> odata = out()->move_block(a1, x5, x4, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x5, x4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x5, x4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma49
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x0, x2, x1)]);
        sort_indices<2,4,5,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: I1375
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), x2.size(), x1.size());
        dgemm_("T", "N", x5.size()*x4.size()*x0.size(), a1.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x4.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x0.size(), a1.size());
  out()->put_block(odata, a1, x5, x4, x0);
}

void Task1132::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1375
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1, x2, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x3.size(), a1.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, a1, x2, x1);
}

void Task1133::Task_local::compute() {
  const Index a1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  // tensor label: I1374
  std::unique_ptr<double[]> odata = out()->move_block(a1, x5, x4, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x5, x4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x5, x4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma240
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x2, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x2, x1, x0)]);
        sort_indices<2,3,4,0,1,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());
        // tensor label: I1381
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, a1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, a1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), a1.size());
        dgemm_("T", "N", x5.size()*x4.size()*x0.size(), a1.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x4.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x0.size(), a1.size());
  out()->put_block(odata, a1, x5, x4, x0);
}

void Task1134::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  // tensor label: I1381
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, a1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, a1);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x3.size(), x2.size(), x1.size(), a1.size());
  }
  out()->put_block(odata, x3, x2, x1, a1);
}

void Task1135::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c2, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c2, x4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c2.size(), x4.size());
      // tensor label: I1377
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x5, x0, x4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x5, x0, x4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, a3.size(), x5.size(), x0.size(), x4.size());
      dgemm_("T", "N", a1.size()*c2.size(), a3.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<2,1,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1136::Task_local::compute() {
  const Index a3 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  // tensor label: I1377
  std::unique_ptr<double[]> odata = out()->move_block(a3, x5, x0, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x5, x0, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x5, x0, x4), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma48
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x3, x4, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x3, x4, x2, x1)]);
        sort_indices<2,4,5,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x3.size(), x4.size(), x2.size(), x1.size());
        // tensor label: I1378
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a3, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), a3.size(), x2.size(), x1.size());
        dgemm_("T", "N", x5.size()*x0.size()*x4.size(), a3.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x0.size()*x4.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x0.size(), x4.size(), a3.size());
  out()->put_block(odata, a3, x5, x0, x4);
}

void Task1137::Task_local::compute() {
  const Index x3 = b(0);
  const Index a3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1378
  std::unique_ptr<double[]> odata = out()->move_block(x3, a3, x2, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, x2, x1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x3.size(), a3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, a3, x2, x1);
}

void Task1138::Task_local::compute() {
  const Index a3 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  // tensor label: I1377
  std::unique_ptr<double[]> odata = out()->move_block(a3, x5, x0, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x5, x0, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x5, x0, x4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma230
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x1, x4, x3, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x1, x4, x3, x2)]);
        sort_indices<2,4,5,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x1.size(), x4.size(), x3.size(), x2.size());
        // tensor label: I1384
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, a3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, a3)]);
        sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), a3.size());
        dgemm_("T", "N", x5.size()*x0.size()*x4.size(), a3.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x0.size()*x4.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x0.size(), x4.size(), a3.size());
  out()->put_block(odata, a3, x5, x0, x4);
}

void Task1139::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a3 = b(3);
  // tensor label: I1384
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, a3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, a3);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x3.size(), x2.size(), x1.size(), a3.size());
  }
  out()->put_block(odata, x3, x2, x1, a3);
}

void Task1140::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c4, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, c4, c2, a1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), c4.size(), c2.size(), a1.size());
      // tensor label: I1386
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c4, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c4, x0, x1)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, a3.size(), c4.size(), x0.size(), x1.size());
      dgemm_("T", "N", c2.size()*a1.size(), a3.size()*x0.size(), c4.size()*x1.size(),
             1.0, i0data_sorted, c4.size()*x1.size(), i1data_sorted, c4.size()*x1.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1141::Task_local::compute() {
  const Index a3 = b(0);
  const Index c4 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1386
  std::unique_ptr<double[]> odata = out()->move_block(a3, c4, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c4, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c4, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma27
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      // tensor label: I1387
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a3, c4, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a3, c4, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a3.size(), c4.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), a3.size()*c4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a3.size(), c4.size());
  out()->put_block(odata, a3, c4, x0, x1);
}

void Task1142::Task_local::compute() {
  const Index x3 = b(0);
  const Index a3 = b(1);
  const Index c4 = b(2);
  const Index x2 = b(3);
  // tensor label: I1387
  std::unique_ptr<double[]> odata = out()->move_block(x3, a3, c4, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c4, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), a3.size(), c4.size(), x2.size());
  }
  out()->put_block(odata, x3, a3, c4, x2);
}

void Task1143::Task_local::compute() {
  const Index a3 = b(0);
  const Index c4 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1386
  std::unique_ptr<double[]> odata = out()->move_block(a3, c4, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c4, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c4, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma29
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I1417
      std::unique_ptr<double[]> i1data = in(1)->get_block(c4, a3, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c4, a3, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c4.size(), a3.size(), x3.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), c4.size()*a3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c4.size(), a3.size());
  out()->put_block(odata, a3, c4, x0, x1);
}

void Task1144::Task_local::compute() {
  const Index c4 = b(0);
  const Index a3 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1417
  std::unique_ptr<double[]> odata = out()->move_block(c4, a3, x3, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a3, x3, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c4.size(), a3.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, c4, a3, x3, x2);
}

void Task1145::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c4, c2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, c4, c2, a3)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), c4.size(), c2.size(), a3.size());
      // tensor label: I1389
      std::unique_ptr<double[]> i1data = in(1)->get_block(a1, c4, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, c4, x0, x1)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, a1.size(), c4.size(), x0.size(), x1.size());
      dgemm_("T", "N", c2.size()*a3.size(), a1.size()*x0.size(), c4.size()*x1.size(),
             1.0, i0data_sorted, c4.size()*x1.size(), i1data_sorted, c4.size()*x1.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), a1.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1146::Task_local::compute() {
  const Index a1 = b(0);
  const Index c4 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1389
  std::unique_ptr<double[]> odata = out()->move_block(a1, c4, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, c4, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, c4, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma27
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      // tensor label: I1390
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, c4, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, c4, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), c4.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), a1.size()*c4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), c4.size());
  out()->put_block(odata, a1, c4, x0, x1);
}

void Task1147::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index c4 = b(2);
  const Index x2 = b(3);
  // tensor label: I1390
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1, c4, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c4, x2);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), a1.size(), c4.size(), x2.size());
  }
  out()->put_block(odata, x3, a1, c4, x2);
}

void Task1148::Task_local::compute() {
  const Index a1 = b(0);
  const Index c4 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1389
  std::unique_ptr<double[]> odata = out()->move_block(a1, c4, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, c4, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, c4, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma29
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I1420
      std::unique_ptr<double[]> i1data = in(1)->get_block(c4, a1, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c4, a1, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c4.size(), a1.size(), x3.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), c4.size()*a1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c4.size(), a1.size());
  out()->put_block(odata, a1, c4, x0, x1);
}

void Task1149::Task_local::compute() {
  const Index c4 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1420
  std::unique_ptr<double[]> odata = out()->move_block(c4, a1, x3, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a1, x3, x2);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c4.size(), a1.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, c4, a1, x3, x2);
}

#endif
