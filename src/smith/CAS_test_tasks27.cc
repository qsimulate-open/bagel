//
// BAGEL - Parallel electron correlation program.
// Filename: CAS_test_tasks27.cc
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


#include <src/smith/CAS_test_tasks27.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CAS_test;

void Task1300::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a2 = b(4);
  // tensor label: I1425
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x2, x1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x2, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x2, x1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, x4, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, x4, x3)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), x4.size(), x3.size());
        // tensor label: I1426
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x0, x4, x3, x2, x1)]);
        sort_indices<4,3,1,0,2,5,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        dgemm_("T", "N", a2.size(), ci0.size()*x0.size()*x2.size()*x1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a2.size(), ci0.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x0, x2, x1, a2);
}

void Task1301::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);
  // tensor label: I1426
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x3, x2, x1);
  {
    // tensor label: Gamma437
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,6,1,1,1,4>(i0data, odata, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x0, x4, x3, x2, x1);
}

void Task1302::Task_local::compute() {
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
    for (auto& a3 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(a3, c2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, c2)]);
      sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a3.size(), c2.size());
      // tensor label: I1429
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x2, x1, a3, c2, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x2, x1, a3, c2, a1)]);
      sort_indices<5,4,0,1,2,3,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x2.size(), x1.size(), a3.size(), c2.size(), a1.size());
      dgemm_("T", "N", 1, ci0.size()*x0.size()*x2.size()*x1.size()*a1.size(), a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,3,4,1,1,1,1>(odata_sorted, odata, ci0.size(), x0.size(), x2.size(), x1.size(), a1.size());
  out()->put_block(odata, ci0, x0, x2, x1, a1);
}

void Task1303::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a3 = b(4);
  const Index c2 = b(5);
  const Index a1 = b(6);
  // tensor label: I1429
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x2, x1, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x2, x1, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x2, x1, a3, c2, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a1.size());
    // tensor label: I1430
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
    sort_indices<1,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), ci0.size()*x0.size()*x2.size()*x1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }
  sort_indices<3,4,5,6,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), ci0.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x0, x2, x1, a3, c2, a1);
}

void Task1304::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: I1430
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma438
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,-1,4>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}

void Task1305::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a3 = b(4);
  const Index c2 = b(5);
  const Index a1 = b(6);
  // tensor label: I1429
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x2, x1, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x2, x1, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x2, x1, a3, c2, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());
    // tensor label: I1434
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
    sort_indices<1,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), ci0.size()*x0.size()*x2.size()*x1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }
  sort_indices<3,4,5,6,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), ci0.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x0, x2, x1, a3, c2, a1);
}

void Task1306::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: I1434
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma438
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,1,2>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}

void Task1307::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a1 = b(4);
  // tensor label: I1413
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x2, x1, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(a2, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a2, x3)]);
      sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a2.size(), x3.size());
      // tensor label: I1437
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x3, x2, x1, a1, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x3, x2, x1, a1, a2)]);
      sort_indices<2,6,0,1,3,4,5,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x3.size(), x2.size(), x1.size(), a1.size(), a2.size());
      dgemm_("T", "N", 1, ci0.size()*x0.size()*x2.size()*x1.size()*a1.size(), x3.size()*a2.size(),
             1.0, i0data_sorted, x3.size()*a2.size(), i1data_sorted, x3.size()*a2.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,3,4,1,1,1,1>(odata_sorted, odata, ci0.size(), x0.size(), x2.size(), x1.size(), a1.size());
  out()->put_block(odata, ci0, x0, x2, x1, a1);
}

void Task1308::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  const Index a1 = b(5);
  const Index a2 = b(6);
  // tensor label: I1437
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x3, x2, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x3, x2, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x3, x2, x1, a1, a2), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a2.size());
      // tensor label: I1438
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x0, x4, x3, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x0, x4, x3, x2, x1)]);
      sort_indices<3,1,0,2,4,5,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
      dgemm_("T", "N", a1.size()*a2.size(), ci0.size()*x0.size()*x3.size()*x2.size()*x1.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }
  sort_indices<2,3,4,5,6,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), ci0.size(), x0.size(), x3.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x0, x3, x2, x1, a1, a2);
}

void Task1309::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);
  // tensor label: I1438
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x3, x2, x1);
  {
    // tensor label: Gamma437
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,6,1,1,1,2>(i0data, odata, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x0, x4, x3, x2, x1);
}

void Task1310::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a1 = b(4);
  // tensor label: I1413
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x2, x1, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, x3)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
        // tensor label: I1946
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x0, x4, x3, x2, x1)]);
        sort_indices<4,3,1,0,2,5,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        dgemm_("T", "N", a1.size(), ci0.size()*x0.size()*x2.size()*x1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x0, x2, x1, a1);
}

void Task1311::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);
  // tensor label: I1946
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x3, x2, x1);
  {
    // tensor label: Gamma437
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x3, x2, x1);
    dscal_(ci0.size()*x5.size()*x0.size()*x4.size()*x3.size()*x2.size()*x1.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,4,5,6,1,1,-1,4>(i0data, odata, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x0, x4, x3, x2, x1);
}

void Task1312::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a1 = b(4);
  // tensor label: I1413
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x2, x1, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, x3)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
        // tensor label: I2033
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x0, x4, x3, x2, x1)]);
        sort_indices<4,3,1,0,2,5,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        dgemm_("T", "N", a1.size(), ci0.size()*x0.size()*x2.size()*x1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x0, x2, x1, a1);
}

void Task1313::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);
  // tensor label: I2033
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x3, x2, x1);
  {
    // tensor label: Gamma437
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,6,1,1,1,2>(i0data, odata, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x0, x4, x3, x2, x1);
}

void Task1314::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a1 = b(4);
  // tensor label: I1413
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x2, x1, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, a1)]);
        sort_indices<2,1,0,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), a1.size());
        // tensor label: I2036
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x4, x3, x0, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x4, x3, x0, x2, x1)]);
        sort_indices<3,2,1,0,4,5,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
        dgemm_("T", "N", a1.size(), ci0.size()*x0.size()*x2.size()*x1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x0, x2, x1, a1);
}

void Task1315::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x0 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);
  // tensor label: I2036
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x3, x0, x2, x1);
  {
    // tensor label: Gamma435
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x4, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,5,6,1,1,1,2>(i0data, odata, ci0.size(), x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x4, x3, x0, x2, x1);
}

void Task1316::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a1 = b(4);
  // tensor label: I1413
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x2, x1, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size());
    // tensor label: I2111
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
    sort_indices<1,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
    dgemm_("T", "N", a1.size(), ci0.size()*x0.size()*x2.size()*x1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x0, x2, x1, a1);
}

void Task1317::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: I2111
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma438
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}

void Task1318::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        for (auto& a4 : *range_[2]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
          // tensor label: I1440
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, c1, a4, c3, a2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, c1, a4, c3, a2)]);
          sort_indices<1,4,3,2,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), c1.size(), a4.size(), c3.size(), a2.size());
          dgemm_("T", "N", 1, ci0.size(), c1.size()*a4.size()*c3.size()*a2.size(),
                 1.0, i0data_sorted, c1.size()*a4.size()*c3.size()*a2.size(), i1data_sorted, c1.size()*a4.size()*c3.size()*a2.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1319::Task_local::compute() {
  const Index ci0 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);
  // tensor label: I1440
  std::unique_ptr<double[]> odata = out()->move_block(ci0, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, c1, a4, c3, a2), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), a2.size());
    // tensor label: I1441
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, c1, a4, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, c1, a4, c3)]);
    sort_indices<1,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), c1.size(), a4.size(), c3.size());
    dgemm_("T", "N", a2.size(), ci0.size()*c1.size()*a4.size()*c3.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a2.size(), ci0.size(), c1.size(), a4.size(), c3.size());
  out()->put_block(odata, ci0, c1, a4, c3, a2);
}

void Task1320::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a4 = b(3);
  const Index c3 = b(4);
  // tensor label: I1441
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c1, a4, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x1.size());
    // tensor label: I1442
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size());
    dgemm_("T", "N", c1.size()*a4.size()*c3.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size()*a4.size()*c3.size());
  }
  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, c1, a4, c3);
}

void Task1321::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  // tensor label: I1442
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1);
  {
    // tensor label: Gamma394
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
    sort_indices<0,1,2,1,1,-1,2>(i0data, odata, ci0.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x0, x1);
}

void Task1322::Task_local::compute() {
  const Index ci0 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);
  // tensor label: I1440
  std::unique_ptr<double[]> odata = out()->move_block(ci0, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, c1, a4, c3, a2), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a4)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), a4.size());
    // tensor label: I1445
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, c1, a2, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, c1, a2, c3)]);
    sort_indices<1,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), c1.size(), a2.size(), c3.size());
    dgemm_("T", "N", a4.size(), ci0.size()*c1.size()*a2.size()*c3.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a4.size());
  }
  sort_indices<1,2,0,4,3,1,1,1,1>(odata_sorted, odata, a4.size(), ci0.size(), c1.size(), a2.size(), c3.size());
  out()->put_block(odata, ci0, c1, a4, c3, a2);
}

void Task1323::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  const Index c3 = b(4);
  // tensor label: I1445
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c1, a2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c1, a2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c1, a2, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x1.size());
    // tensor label: I1446
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size());
    dgemm_("T", "N", c1.size()*a2.size()*c3.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, c1, a2, c3);
}

void Task1324::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  // tensor label: I1446
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1);
  {
    // tensor label: Gamma394
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
    sort_indices<0,1,2,1,1,1,1>(i0data, odata, ci0.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x0, x1);
}

void Task1325::Task_local::compute() {
  const Index ci0 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);
  // tensor label: I1440
  std::unique_ptr<double[]> odata = out()->move_block(ci0, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, c1, a4, c3, a2), 0.0);
  // tensor label: f1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size());
  // tensor label: I1449
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, a4, c1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, a4, c1)]);
  sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), a4.size(), c1.size());
  dgemm_("T", "N", c3.size()*a2.size(), ci0.size()*a4.size()*c1.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c3.size()*a2.size());
  sort_indices<2,4,3,0,1,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), ci0.size(), a4.size(), c1.size());
  out()->put_block(odata, ci0, c1, a4, c3, a2);
}

void Task1326::Task_local::compute() {
  const Index ci0 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  // tensor label: I1449
  std::unique_ptr<double[]> odata = out()->move_block(ci0, a4, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, a4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, a4, c1), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c1, x0)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c1.size(), x0.size());
      // tensor label: I1450
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
      sort_indices<2,1,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
      dgemm_("T", "N", a4.size()*c1.size(), ci0.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, a4.size()*c1.size());
    }
  }
  sort_indices<2,0,1,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), ci0.size());
  out()->put_block(odata, ci0, a4, c1);
}

void Task1327::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1450
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1328::Task_local::compute() {
  const Index ci0 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  // tensor label: I1449
  std::unique_ptr<double[]> odata = out()->move_block(ci0, a4, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, a4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, a4, c1), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, x1, x0)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), x1.size(), x0.size());
      // tensor label: I1458
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
      sort_indices<2,1,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a4.size(), ci0.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, c1.size()*a4.size());
    }
  }
  sort_indices<2,1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), ci0.size());
  out()->put_block(odata, ci0, a4, c1);
}

void Task1329::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1458
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1330::Task_local::compute() {
  const Index ci0 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  // tensor label: I1449
  std::unique_ptr<double[]> odata = out()->move_block(ci0, a4, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, a4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, a4, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a4, c1, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a4, c1, x1)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), a4.size(), c1.size(), x1.size());
      // tensor label: I1812
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
      sort_indices<1,2,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
      dgemm_("T", "N", a4.size()*c1.size(), ci0.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, a4.size()*c1.size());
    }
  }
  sort_indices<2,0,1,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), ci0.size());
  out()->put_block(odata, ci0, a4, c1);
}

void Task1331::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1812
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1332::Task_local::compute() {
  const Index ci0 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  // tensor label: I1449
  std::unique_ptr<double[]> odata = out()->move_block(ci0, a4, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, a4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, a4, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, x0, x1)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), x0.size(), x1.size());
      // tensor label: I1820
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
      sort_indices<1,2,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a4.size(), ci0.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, c1.size()*a4.size());
    }
  }
  sort_indices<2,1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), ci0.size());
  out()->put_block(odata, ci0, a4, c1);
}

void Task1333::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1820
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1334::Task_local::compute() {
  const Index ci0 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);
  // tensor label: I1440
  std::unique_ptr<double[]> odata = out()->move_block(ci0, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, c1, a4, c3, a2), 0.0);
  // tensor label: f1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size());
  // tensor label: I1453
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, a2, c1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, a2, c1)]);
  sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), a2.size(), c1.size());
  dgemm_("T", "N", c3.size()*a4.size(), ci0.size()*a2.size()*c1.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c3.size()*a4.size());
  sort_indices<2,4,1,0,3,1,1,1,1>(odata_sorted, odata, c3.size(), a4.size(), ci0.size(), a2.size(), c1.size());
  out()->put_block(odata, ci0, c1, a4, c3, a2);
}

void Task1335::Task_local::compute() {
  const Index ci0 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  // tensor label: I1453
  std::unique_ptr<double[]> odata = out()->move_block(ci0, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, a2, c1), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2, c1, x0)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), x0.size());
      // tensor label: I1454
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
      sort_indices<2,1,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
      dgemm_("T", "N", a2.size()*c1.size(), ci0.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<2,0,1,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), ci0.size());
  out()->put_block(odata, ci0, a2, c1);
}

void Task1336::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1454
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1337::Task_local::compute() {
  const Index ci0 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  // tensor label: I1453
  std::unique_ptr<double[]> odata = out()->move_block(ci0, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, a2, c1), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x1, x0)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x1.size(), x0.size());
      // tensor label: I1462
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
      sort_indices<2,1,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), ci0.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<2,1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), ci0.size());
  out()->put_block(odata, ci0, a2, c1);
}

void Task1338::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1462
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,2,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1339::Task_local::compute() {
  const Index ci0 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  // tensor label: I1453
  std::unique_ptr<double[]> odata = out()->move_block(ci0, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, a2, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a2, c1, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a2, c1, x1)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), a2.size(), c1.size(), x1.size());
      // tensor label: I1816
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
      sort_indices<1,2,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
      dgemm_("T", "N", a2.size()*c1.size(), ci0.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<2,0,1,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), ci0.size());
  out()->put_block(odata, ci0, a2, c1);
}

void Task1340::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1816
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1341::Task_local::compute() {
  const Index ci0 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  // tensor label: I1453
  std::unique_ptr<double[]> odata = out()->move_block(ci0, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, a2, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x0, x1)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x0.size(), x1.size());
      // tensor label: I1824
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
      sort_indices<1,2,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), ci0.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<2,1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), ci0.size());
  out()->put_block(odata, ci0, a2, c1);
}

void Task1342::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1824
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,2,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1343::Task_local::compute() {
  const Index ci0 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);
  // tensor label: I1440
  std::unique_ptr<double[]> odata = out()->move_block(ci0, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, c1, a4, c3, a2), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
  // tensor label: I1465
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0)]);
  sort_indices<0,0,1,1,1>(i1data, i1data_sorted, ci0.size());
  dgemm_("T", "N", c1.size()*a4.size()*c3.size()*a2.size(), ci0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a4.size()*c3.size()*a2.size());
  sort_indices<4,0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size(), ci0.size());
  out()->put_block(odata, ci0, c1, a4, c3, a2);
}

void Task1344::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1465
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  {
    // tensor label: Gamma447
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0);
    sort_indices<0,1,1,-1,1>(i0data, odata, ci0.size());
  }
  {
    // tensor label: Gamma447
    std::unique_ptr<double[]> i1data = in(0)->get_block(ci0);
    sort_indices<0,1,1,-1,1>(i1data, odata, ci0.size());
  }
  out()->put_block(odata, ci0);
}

void Task1345::Task_local::compute() {
  const Index ci0 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);
  // tensor label: I1440
  std::unique_ptr<double[]> odata = out()->move_block(ci0, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, c1, a4, c3, a2), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
  // tensor label: I1468
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0)]);
  sort_indices<0,0,1,1,1>(i1data, i1data_sorted, ci0.size());
  dgemm_("T", "N", c1.size()*a2.size()*c3.size()*a4.size(), ci0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a2.size()*c3.size()*a4.size());
  sort_indices<4,0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size(), ci0.size());
  out()->put_block(odata, ci0, c1, a4, c3, a2);
}

void Task1346::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1468
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  {
    // tensor label: Gamma447
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0);
    sort_indices<0,1,1,2,1>(i0data, odata, ci0.size());
  }
  {
    // tensor label: Gamma447
    std::unique_ptr<double[]> i1data = in(0)->get_block(ci0);
    sort_indices<0,1,1,2,1>(i1data, odata, ci0.size());
  }
  out()->put_block(odata, ci0);
}

void Task1347::Task_local::compute() {
  const Index ci0 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);
  // tensor label: I1440
  std::unique_ptr<double[]> odata = out()->move_block(ci0, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, c1, a4, c3, a2), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), x0.size());
    // tensor label: I1471
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, a4, c1, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, a4, c1, a2)]);
    sort_indices<1,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), a4.size(), c1.size(), a2.size());
    dgemm_("T", "N", c3.size(), ci0.size()*a4.size()*c1.size()*a2.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, c3.size());
  }
  sort_indices<1,3,2,0,4,1,1,1,1>(odata_sorted, odata, c3.size(), ci0.size(), a4.size(), c1.size(), a2.size());
  out()->put_block(odata, ci0, c1, a4, c3, a2);
}

void Task1348::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a4 = b(2);
  const Index c1 = b(3);
  const Index a2 = b(4);
  // tensor label: I1471
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a4, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a4, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a4, c1, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c1, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c1.size(), a2.size());
    // tensor label: I1472
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
    dgemm_("T", "N", a4.size()*c1.size()*a2.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a4.size()*c1.size()*a2.size());
  }
  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), a2.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, a4, c1, a2);
}

void Task1349::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1472
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

