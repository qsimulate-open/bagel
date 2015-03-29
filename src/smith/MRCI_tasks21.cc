//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks21.cc
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

#include <src/smith/MRCI_tasks21.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task1000::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& c5 : *range_[0]) {
    for (auto& a6 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c5, a6);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c5, a6)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c5.size(), a6.size());
      // tensor label: I1304
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, a6, c5);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, a6, c5)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), a6.size(), c5.size());
      dgemm_("T", "N", c1.size()*a4.size(), c3.size()*a2.size(), a6.size()*c5.size(),
             1.0, i0data_sorted, a6.size()*c5.size(), i1data_sorted, a6.size()*c5.size(),
             1.0, odata_sorted, c1.size()*a4.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task1001::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index a6 = b(2);
  const Index c5 = b(3);
  // tensor label: I1304
  std::unique_ptr<double[]> odata = out()->move_block(c3, a2, a6, c5);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, a6, c5);
    sort_indices<0,1,2,3,1,1,-8,1>(i0data, odata, c3.size(), a2.size(), a6.size(), c5.size());
  }
  out()->put_block(odata, c3, a2, a6, c5);
}

void Task1002::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& a6 : *range_[2]) {
    for (auto& c5 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a6, c5, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a6, c5, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a6.size(), c5.size(), a2.size());
      // tensor label: I1306
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a4, a6, c5);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a4, a6, c5)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), a4.size(), a6.size(), c5.size());
      dgemm_("T", "N", c1.size()*a2.size(), c3.size()*a4.size(), a6.size()*c5.size(),
             1.0, i0data_sorted, a6.size()*c5.size(), i1data_sorted, a6.size()*c5.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task1003::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index a6 = b(2);
  const Index c5 = b(3);
  // tensor label: I1306
  std::unique_ptr<double[]> odata = out()->move_block(c3, a4, a6, c5);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, a6, c5);
    sort_indices<0,1,2,3,1,1,-8,1>(i0data, odata, c3.size(), a4.size(), a6.size(), c5.size());
  }
  out()->put_block(odata, c3, a4, a6, c5);
}

void Task1004::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& c5 : *range_[0]) {
    for (auto& a6 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c5, a6);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c5, a6)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c5.size(), a6.size());
      // tensor label: I1308
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a4, a6, c5);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a4, a6, c5)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), a4.size(), a6.size(), c5.size());
      dgemm_("T", "N", c1.size()*a2.size(), c3.size()*a4.size(), a6.size()*c5.size(),
             1.0, i0data_sorted, a6.size()*c5.size(), i1data_sorted, a6.size()*c5.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task1005::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index a6 = b(2);
  const Index c5 = b(3);
  // tensor label: I1308
  std::unique_ptr<double[]> odata = out()->move_block(c3, a4, a6, c5);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, a6, c5);
    sort_indices<0,1,2,3,1,1,16,1>(i0data, odata, c3.size(), a4.size(), a6.size(), c5.size());
  }
  out()->put_block(odata, c3, a4, a6, c5);
}

void Task1006::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a4, c1, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c1.size(), a2.size());
    // tensor label: I1314
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), x3.size());
    dgemm_("T", "N", a4.size()*c1.size()*a2.size(), c3.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a4.size()*c1.size()*a2.size());
  }
  sort_indices<2,1,0,3,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), a2.size(), c3.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task1007::Task_local::compute() {
  const Index c3 = b(0);
  const Index x3 = b(1);
  // tensor label: I1314
  std::unique_ptr<double[]> odata = out()->move_block(c3, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: Gamma29
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
        // tensor label: I1315
        std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x2, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x2, x1, x0)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), x2.size(), x1.size(), x0.size());
        dgemm_("T", "N", x3.size(), c3.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), c3.size());
  out()->put_block(odata, c3, x3);
}

void Task1008::Task_local::compute() {
  const Index c3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1315
  std::unique_ptr<double[]> odata = out()->move_block(c3, x2, x1, x0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, c3, x2, x1, x0);
}

void Task1009::Task_local::compute() {
  const Index c3 = b(0);
  const Index x3 = b(1);
  // tensor label: I1314
  std::unique_ptr<double[]> odata = out()->move_block(c3, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x3), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma51
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: I1321
        std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, c3, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, c3, x0)]);
        sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), c3.size(), x0.size());
        dgemm_("T", "N", x3.size(), c3.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), c3.size());
  out()->put_block(odata, c3, x3);
}

void Task1010::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index c3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1321
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, c3, x0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x1, c3, x0);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x2.size(), x1.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x2, x1, c3, x0);
}

void Task1011::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c1, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c1.size(), a4.size());
    // tensor label: I1317
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), x3.size());
    dgemm_("T", "N", a2.size()*c1.size()*a4.size(), c3.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a2.size()*c1.size()*a4.size());
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), a4.size(), c3.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task1012::Task_local::compute() {
  const Index c3 = b(0);
  const Index x3 = b(1);
  // tensor label: I1317
  std::unique_ptr<double[]> odata = out()->move_block(c3, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: Gamma29
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
        // tensor label: I1318
        std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x2, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x2, x1, x0)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), x2.size(), x1.size(), x0.size());
        dgemm_("T", "N", x3.size(), c3.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), c3.size());
  out()->put_block(odata, c3, x3);
}

void Task1013::Task_local::compute() {
  const Index c3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1318
  std::unique_ptr<double[]> odata = out()->move_block(c3, x2, x1, x0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, c3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, c3, x2, x1, x0);
}

void Task1014::Task_local::compute() {
  const Index c3 = b(0);
  const Index x3 = b(1);
  // tensor label: I1317
  std::unique_ptr<double[]> odata = out()->move_block(c3, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x3), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma51
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: I1324
        std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, c3, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, c3, x0)]);
        sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), c3.size(), x0.size());
        dgemm_("T", "N", x3.size(), c3.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), c3.size());
  out()->put_block(odata, c3, x3);
}

void Task1015::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index c3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1324
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, c3, x0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x1, c3, x0);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x2.size(), x1.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x2, x1, c3, x0);
}

void Task1016::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& c5 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c5, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c5, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c5.size(), a2.size());
      // tensor label: I1326
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c3, c5, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c3, c5, x1)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), c3.size(), c5.size(), x1.size());
      dgemm_("T", "N", a4.size()*a2.size(), c1.size()*c3.size(), c5.size()*x1.size(),
             1.0, i0data_sorted, c5.size()*x1.size(), i1data_sorted, c5.size()*x1.size(),
             1.0, odata_sorted, a4.size()*a2.size());
    }
  }
  sort_indices<1,2,0,3,1,1,1,1>(odata_sorted, odata, a4.size(), a2.size(), c1.size(), c3.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task1017::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index c5 = b(2);
  const Index x1 = b(3);
  // tensor label: I1326
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, c5, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, c5, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, c5, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma32
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I1327
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x0, c3, c5);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x0, c3, c5)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), x0.size(), c3.size(), c5.size());
    dgemm_("T", "N", x1.size(), c1.size()*c3.size()*c5.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), c1.size(), c3.size(), c5.size());
  out()->put_block(odata, c1, c3, c5, x1);
}

void Task1018::Task_local::compute() {
  const Index c1 = b(0);
  const Index x0 = b(1);
  const Index c3 = b(2);
  const Index c5 = b(3);
  // tensor label: I1327
  std::unique_ptr<double[]> odata = out()->move_block(c1, x0, c3, c5);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x0, c3, c5);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c1.size(), x0.size(), c3.size(), c5.size());
  }
  out()->put_block(odata, c1, x0, c3, c5);
}

void Task1019::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& c5 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c5, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2, c5, a4)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c5.size(), a4.size());
      // tensor label: I1329
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c3, c5, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c3, c5, x1)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), c3.size(), c5.size(), x1.size());
      dgemm_("T", "N", a2.size()*a4.size(), c1.size()*c3.size(), c5.size()*x1.size(),
             1.0, i0data_sorted, c5.size()*x1.size(), i1data_sorted, c5.size()*x1.size(),
             1.0, odata_sorted, a2.size()*a4.size());
    }
  }
  sort_indices<0,2,1,3,1,1,1,1>(odata_sorted, odata, a2.size(), a4.size(), c1.size(), c3.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task1020::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index c5 = b(2);
  const Index x1 = b(3);
  // tensor label: I1329
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, c5, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, c5, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, c5, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma32
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I1330
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x0, c3, c5);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x0, c3, c5)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), x0.size(), c3.size(), c5.size());
    dgemm_("T", "N", x1.size(), c1.size()*c3.size()*c5.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), c1.size(), c3.size(), c5.size());
  out()->put_block(odata, c1, c3, c5, x1);
}

void Task1021::Task_local::compute() {
  const Index c1 = b(0);
  const Index x0 = b(1);
  const Index c3 = b(2);
  const Index c5 = b(3);
  // tensor label: I1330
  std::unique_ptr<double[]> odata = out()->move_block(c1, x0, c3, c5);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x0, c3, c5);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, c1.size(), x0.size(), c3.size(), c5.size());
  }
  out()->put_block(odata, c1, x0, c3, c5);
}

void Task1022::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a5, c1, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a5, c1, a4)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a5.size(), c1.size(), a4.size());
      // tensor label: I1332
      std::unique_ptr<double[]> i1data = in(1)->get_block(a5, c3, a2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, c3, a2, x1)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a5.size(), c3.size(), a2.size(), x1.size());
      dgemm_("T", "N", c1.size()*a4.size(), c3.size()*a2.size(), a5.size()*x1.size(),
             1.0, i0data_sorted, a5.size()*x1.size(), i1data_sorted, a5.size()*x1.size(),
             1.0, odata_sorted, c1.size()*a4.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task1023::Task_local::compute() {
  const Index a5 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1332
  std::unique_ptr<double[]> odata = out()->move_block(a5, c3, a2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a5, c3, a2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a5, c3, a2, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma32
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I1333
    std::unique_ptr<double[]> i1data = in(1)->get_block(a5, x0, c3, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, x0, c3, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a5.size(), x0.size(), c3.size(), a2.size());
    dgemm_("T", "N", x1.size(), a5.size()*c3.size()*a2.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), a5.size(), c3.size(), a2.size());
  out()->put_block(odata, a5, c3, a2, x1);
}

void Task1024::Task_local::compute() {
  const Index a5 = b(0);
  const Index x0 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I1333
  std::unique_ptr<double[]> odata = out()->move_block(a5, x0, c3, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, x0, c3, a2);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, a5.size(), x0.size(), c3.size(), a2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c3, x0, a5, a2);
    sort_indices<2,1,0,3,1,1,1,1>(i1data, odata, c3.size(), x0.size(), a5.size(), a2.size());
  }
  out()->put_block(odata, a5, x0, c3, a2);
}

void Task1025::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c1, a5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c1, a5)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c1.size(), a5.size());
      // tensor label: I1335
      std::unique_ptr<double[]> i1data = in(1)->get_block(a5, c3, a2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, c3, a2, x1)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a5.size(), c3.size(), a2.size(), x1.size());
      dgemm_("T", "N", a4.size()*c1.size(), c3.size()*a2.size(), a5.size()*x1.size(),
             1.0, i0data_sorted, a5.size()*x1.size(), i1data_sorted, a5.size()*x1.size(),
             1.0, odata_sorted, a4.size()*c1.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), c3.size(), a2.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task1026::Task_local::compute() {
  const Index a5 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1335
  std::unique_ptr<double[]> odata = out()->move_block(a5, c3, a2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a5, c3, a2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a5, c3, a2, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma32
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I1336
    std::unique_ptr<double[]> i1data = in(1)->get_block(a5, x0, c3, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, x0, c3, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a5.size(), x0.size(), c3.size(), a2.size());
    dgemm_("T", "N", x1.size(), a5.size()*c3.size()*a2.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), a5.size(), c3.size(), a2.size());
  out()->put_block(odata, a5, c3, a2, x1);
}

void Task1027::Task_local::compute() {
  const Index a5 = b(0);
  const Index x0 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I1336
  std::unique_ptr<double[]> odata = out()->move_block(a5, x0, c3, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, x0, c3, a2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, a5.size(), x0.size(), c3.size(), a2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c3, x0, a5, a2);
    sort_indices<2,1,0,3,1,1,-2,1>(i1data, odata, c3.size(), x0.size(), a5.size(), a2.size());
  }
  out()->put_block(odata, a5, x0, c3, a2);
}

void Task1028::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a5, c1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a5, c1, a2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a5.size(), c1.size(), a2.size());
      // tensor label: I1338
      std::unique_ptr<double[]> i1data = in(1)->get_block(a5, c3, a4, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, c3, a4, x1)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a5.size(), c3.size(), a4.size(), x1.size());
      dgemm_("T", "N", c1.size()*a2.size(), c3.size()*a4.size(), a5.size()*x1.size(),
             1.0, i0data_sorted, a5.size()*x1.size(), i1data_sorted, a5.size()*x1.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task1029::Task_local::compute() {
  const Index a5 = b(0);
  const Index c3 = b(1);
  const Index a4 = b(2);
  const Index x1 = b(3);
  // tensor label: I1338
  std::unique_ptr<double[]> odata = out()->move_block(a5, c3, a4, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a5, c3, a4, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a5, c3, a4, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma32
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I1339
    std::unique_ptr<double[]> i1data = in(1)->get_block(a5, x0, c3, a4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, x0, c3, a4)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a5.size(), x0.size(), c3.size(), a4.size());
    dgemm_("T", "N", x1.size(), a5.size()*c3.size()*a4.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), a5.size(), c3.size(), a4.size());
  out()->put_block(odata, a5, c3, a4, x1);
}

void Task1030::Task_local::compute() {
  const Index a5 = b(0);
  const Index x0 = b(1);
  const Index c3 = b(2);
  const Index a4 = b(3);
  // tensor label: I1339
  std::unique_ptr<double[]> odata = out()->move_block(a5, x0, c3, a4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, x0, c3, a4);
    sort_indices<0,1,2,3,1,1,4,1>(i0data, odata, a5.size(), x0.size(), c3.size(), a4.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c3, x0, a5, a4);
    sort_indices<2,1,0,3,1,1,-2,1>(i1data, odata, c3.size(), x0.size(), a5.size(), a4.size());
  }
  out()->put_block(odata, a5, x0, c3, a4);
}

void Task1031::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c1, a5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2, c1, a5)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), a5.size());
      // tensor label: I1341
      std::unique_ptr<double[]> i1data = in(1)->get_block(a5, c3, a4, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, c3, a4, x1)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a5.size(), c3.size(), a4.size(), x1.size());
      dgemm_("T", "N", a2.size()*c1.size(), c3.size()*a4.size(), a5.size()*x1.size(),
             1.0, i0data_sorted, a5.size()*x1.size(), i1data_sorted, a5.size()*x1.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), c3.size(), a4.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task1032::Task_local::compute() {
  const Index a5 = b(0);
  const Index c3 = b(1);
  const Index a4 = b(2);
  const Index x1 = b(3);
  // tensor label: I1341
  std::unique_ptr<double[]> odata = out()->move_block(a5, c3, a4, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a5, c3, a4, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a5, c3, a4, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma32
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I1342
    std::unique_ptr<double[]> i1data = in(1)->get_block(a5, x0, c3, a4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, x0, c3, a4)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a5.size(), x0.size(), c3.size(), a4.size());
    dgemm_("T", "N", x1.size(), a5.size()*c3.size()*a4.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), a5.size(), c3.size(), a4.size());
  out()->put_block(odata, a5, c3, a4, x1);
}

void Task1033::Task_local::compute() {
  const Index a5 = b(0);
  const Index x0 = b(1);
  const Index c3 = b(2);
  const Index a4 = b(3);
  // tensor label: I1342
  std::unique_ptr<double[]> odata = out()->move_block(a5, x0, c3, a4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, x0, c3, a4);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, a5.size(), x0.size(), c3.size(), a4.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c3, x0, a5, a4);
    sort_indices<2,1,0,3,1,1,1,1>(i1data, odata, c3.size(), x0.size(), a5.size(), a4.size());
  }
  out()->put_block(odata, a5, x0, c3, a4);
}

void Task1034::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c2, a3, x0, a1);
  {
    // tensor label: I194
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, c2, x0, a1);
    sort_indices<1,0,2,3,1,1,1,1>(i0data, odata, a3.size(), c2.size(), x0.size(), a1.size());
  }
  out()->put_block(odata, c2, a3, x0, a1);
}

void Task1035::Task_local::compute() {
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
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size());
    // tensor label: I195
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, x1, x0)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), x1.size(), x0.size());
    dgemm_("T", "N", a1.size(), a3.size()*c2.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), c2.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1036::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I195
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma29
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I196
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a3, c2, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a3, c2, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a3.size(), c2.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), a3.size()*c2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a3.size(), c2.size());
  out()->put_block(odata, a3, c2, x1, x0);
}

void Task1037::Task_local::compute() {
  const Index x3 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index x2 = b(3);
  // tensor label: I196
  std::unique_ptr<double[]> odata = out()->move_block(x3, a3, c2, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), a3.size(), c2.size(), x2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c2, a3, x3, x2);
    sort_indices<2,1,0,3,1,1,2,1>(i1data, odata, c2.size(), a3.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, x3, a3, c2, x2);
}

void Task1038::Task_local::compute() {
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
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size());
    // tensor label: I198
    std::unique_ptr<double[]> i1data = in(1)->get_block(a1, c2, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, c2, x0, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a1.size(), c2.size(), x0.size(), x1.size());
    dgemm_("T", "N", a3.size(), a1.size()*c2.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size());
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, a3.size(), a1.size(), c2.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1039::Task_local::compute() {
  const Index a1 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I198
  std::unique_ptr<double[]> odata = out()->move_block(a1, c2, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, c2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, c2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma27
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      // tensor label: I199
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, c2, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, c2, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), c2.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), a1.size()*c2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), c2.size());
  out()->put_block(odata, a1, c2, x0, x1);
}

void Task1040::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index c2 = b(2);
  const Index x2 = b(3);
  // tensor label: I199
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1, c2, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), a1.size(), c2.size(), x2.size());
  }
  out()->put_block(odata, x3, a1, c2, x2);
}

void Task1041::Task_local::compute() {
  const Index a1 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I198
  std::unique_ptr<double[]> odata = out()->move_block(a1, c2, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, c2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, c2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma29
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I205
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x3.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), c2.size()*a1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a1.size());
  out()->put_block(odata, a1, c2, x0, x1);
}

void Task1042::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I205
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, x3, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x3, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c2.size(), a1.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, c2, a1, x3, x2);
}

void Task1043::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  // tensor label: h1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size());
  // tensor label: I207
  std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), x0.size());
  dgemm_("T", "N", c2.size()*a1.size(), a3.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a1.size());
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1044::Task_local::compute() {
  const Index a3 = b(0);
  const Index x0 = b(1);
  // tensor label: I207
  std::unique_ptr<double[]> odata = out()->move_block(a3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma51
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: I208
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a3, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), a3.size(), x2.size(), x1.size());
        dgemm_("T", "N", x0.size(), a3.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a3.size());
  out()->put_block(odata, a3, x0);
}

void Task1045::Task_local::compute() {
  const Index x3 = b(0);
  const Index a3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I208
  std::unique_ptr<double[]> odata = out()->move_block(x3, a3, x2, x1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), a3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, a3, x2, x1);
}

void Task1046::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  // tensor label: h1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size());
  // tensor label: I210
  std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a1.size(), x0.size());
  dgemm_("T", "N", c2.size()*a3.size(), a1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a3.size());
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), a1.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1047::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  // tensor label: I210
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma51
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: I211
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), x2.size(), x1.size());
        dgemm_("T", "N", x0.size(), a1.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a1.size());
  out()->put_block(odata, a1, x0);
}

void Task1048::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I211
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1, x2, x1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, x1);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), a1.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, a1, x2, x1);
}

void Task1049::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c4, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c4, a1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), a1.size());
    // tensor label: I213
    std::unique_ptr<double[]> i1data = in(1)->get_block(c4, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c4, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c4.size(), x0.size());
    dgemm_("T", "N", c2.size()*a3.size()*a1.size(), x0.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c2.size()*a3.size()*a1.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), a1.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

#endif
