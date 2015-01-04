//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks28.cc
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


#include <src/smith/CASPT2_tasks28.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task1350::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1309
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, x1, x0), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a2, c3, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a2, c3, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a2.size(), c3.size(), x1.size());
      // tensor label: I1644
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c3, a2, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c3, a2, x2)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), c3.size(), a2.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), x3.size()*x2.size(), c3.size()*a2.size(),
             1.0, i0data_sorted, c3.size()*a2.size(), i1data_sorted, c3.size()*a2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x3.size(), x2.size());
  out()->put_block(odata, x2, x3, x1, x0);
}

void Task1351::Task_local::compute() {
  const Index x3 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I1644
  std::unique_ptr<double[]> odata = out()->move_block(x3, c3, a2, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, c3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, c3, a2, x2), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), x2.size());
    // tensor label: I1645
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c3, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c3, a2, c1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), c3.size(), a2.size(), c1.size());
    dgemm_("T", "N", x2.size(), x3.size()*c3.size()*a2.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x2.size(), x3.size(), c3.size(), a2.size());
  out()->put_block(odata, x3, c3, a2, x2);
}

void Task1352::Task_local::compute() {
  const Index x3 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1645
  std::unique_ptr<double[]> odata = out()->move_block(x3, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x3);
    sort_indices<3,2,1,0,1,1,1,4>(i0data, odata, c1.size(), a2.size(), c3.size(), x3.size());
  }
  out()->put_block(odata, x3, c3, a2, c1);
}

void Task1353::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1309
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, x1, x0), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x3, x2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x3, x2, a1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, c2.size(), x3.size(), x2.size(), a1.size());
      // tensor label: I2012
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c2, a1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c2, a1, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), c2.size(), a1.size(), x0.size());
      dgemm_("T", "N", x3.size()*x2.size(), x1.size()*x0.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x3, x1, x0);
}

void Task1354::Task_local::compute() {
  const Index x1 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I2012
  std::unique_ptr<double[]> odata = out()->move_block(x1, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, x0.size(), a1.size(), c2.size(), x1.size());
  }
  out()->put_block(odata, x1, c2, a1, x0);
}

void Task1355::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        for (auto& x4 : *range_[1]) {
          // tensor label: Gamma409
          std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x1, x4);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x5, x0, x1, x4)]);
          sort_indices<1,2,3,4,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x5.size(), x0.size(), x1.size(), x4.size());
          // tensor label: I1317
          std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, x5, x4);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, x5, x4)]);
          sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), x5.size(), x4.size());
          dgemm_("T", "N", ci0.size(), 1, x1.size()*x0.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x1.size()*x0.size()*x5.size()*x4.size(), i1data_sorted, x1.size()*x0.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, ci0.size());
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1356::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  // tensor label: I1317
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x5, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x5, x4), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c2, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c2, x4)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c2.size(), x4.size());
      // tensor label: I1318
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c2, a1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c2, a1, x0)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), c2.size(), a1.size(), x0.size());
      dgemm_("T", "N", x5.size()*x4.size(), x1.size()*x0.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, x5.size()*x4.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, x5, x4);
}

void Task1357::Task_local::compute() {
  const Index x1 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1318
  std::unique_ptr<double[]> odata = out()->move_block(x1, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
    sort_indices<3,2,1,0,1,1,1,4>(i0data, odata, x0.size(), a1.size(), c2.size(), x1.size());
  }
  out()->put_block(odata, x1, c2, a1, x0);
}

void Task1358::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  // tensor label: I1317
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x5, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x5, x4), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), x1.size());
      // tensor label: I1680
      std::unique_ptr<double[]> i1data = in(1)->get_block(x4, c2, a1, x5);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, c2, a1, x5)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x4.size(), c2.size(), a1.size(), x5.size());
      dgemm_("T", "N", x0.size()*x1.size(), x4.size()*x5.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x4.size(), x5.size());
  out()->put_block(odata, x1, x0, x5, x4);
}

void Task1359::Task_local::compute() {
  const Index x4 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x5 = b(3);
  // tensor label: I1680
  std::unique_ptr<double[]> odata = out()->move_block(x4, c2, a1, x5);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c2, x4);
    sort_indices<3,2,1,0,1,1,1,4>(i0data, odata, x5.size(), a1.size(), c2.size(), x4.size());
  }
  out()->put_block(odata, x4, c2, a1, x5);
}

void Task1360::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: Gamma410
          std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x1, x2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x3, x0, x1, x2)]);
          sort_indices<1,2,3,4,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x3.size(), x0.size(), x1.size(), x2.size());
          // tensor label: I1320
          std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, x3, x2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, x3, x2)]);
          sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), x3.size(), x2.size());
          dgemm_("T", "N", ci0.size(), 1, x1.size()*x0.size()*x3.size()*x2.size(),
                 1.0, i0data_sorted, x1.size()*x0.size()*x3.size()*x2.size(), i1data_sorted, x1.size()*x0.size()*x3.size()*x2.size(),
                 1.0, odata_sorted, ci0.size());
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1361::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1320
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c3, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c3.size(), x2.size());
      // tensor label: I1321
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a1, x0, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a1, x0, c3)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), a1.size(), x0.size(), c3.size());
      dgemm_("T", "N", x3.size()*x2.size(), x1.size()*x0.size(), a1.size()*c3.size(),
             1.0, i0data_sorted, a1.size()*c3.size(), i1data_sorted, a1.size()*c3.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, x3, x2);
}

void Task1362::Task_local::compute() {
  const Index x1 = b(0);
  const Index a1 = b(1);
  const Index x0 = b(2);
  const Index c3 = b(3);
  // tensor label: I1321
  std::unique_ptr<double[]> odata = out()->move_block(x1, a1, x0, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, a1, x0, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, a1, x0, c3), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, c3)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), c3.size());
    // tensor label: I1322
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c2, a1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c2, a1, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), c2.size(), a1.size(), x0.size());
    dgemm_("T", "N", c3.size(), x1.size()*a1.size()*x0.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, c3.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c3.size(), x1.size(), a1.size(), x0.size());
  out()->put_block(odata, x1, a1, x0, c3);
}

void Task1363::Task_local::compute() {
  const Index x1 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1322
  std::unique_ptr<double[]> odata = out()->move_block(x1, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
    sort_indices<3,2,1,0,1,1,-1,4>(i0data, odata, x0.size(), a1.size(), c2.size(), x1.size());
  }
  out()->put_block(odata, x1, c2, a1, x0);
}

void Task1364::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1320
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), x2.size());
      // tensor label: I1325
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c2, x0, a3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c2, x0, a3)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), c2.size(), x0.size(), a3.size());
      dgemm_("T", "N", x3.size()*x2.size(), x1.size()*x0.size(), c2.size()*a3.size(),
             1.0, i0data_sorted, c2.size()*a3.size(), i1data_sorted, c2.size()*a3.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, x3, x2);
}

void Task1365::Task_local::compute() {
  const Index x1 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a3 = b(3);
  // tensor label: I1325
  std::unique_ptr<double[]> odata = out()->move_block(x1, c2, x0, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, c2, x0, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, c2, x0, a3), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, a1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a3.size(), a1.size());
    // tensor label: I1326
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c2, a1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c2, a1, x0)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), c2.size(), a1.size(), x0.size());
    dgemm_("T", "N", a3.size(), x1.size()*c2.size()*x0.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, a3.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a3.size(), x1.size(), c2.size(), x0.size());
  out()->put_block(odata, x1, c2, x0, a3);
}

void Task1366::Task_local::compute() {
  const Index x1 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1326
  std::unique_ptr<double[]> odata = out()->move_block(x1, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
    sort_indices<3,2,1,0,1,1,1,4>(i0data, odata, x0.size(), a1.size(), c2.size(), x1.size());
  }
  out()->put_block(odata, x1, c2, a1, x0);
}

void Task1367::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1320
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, x1)]);
      sort_indices<2,1,0,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), x1.size());
      // tensor label: I1356
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x3, a1, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x3, a1, c2)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, x2.size(), x3.size(), a1.size(), c2.size());
      dgemm_("T", "N", x1.size()*x0.size(), x2.size()*x3.size(), a1.size()*c2.size(),
             1.0, i0data_sorted, a1.size()*c2.size(), i1data_sorted, a1.size()*c2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size());
  out()->put_block(odata, x1, x0, x3, x2);
}

void Task1368::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index c2 = b(3);
  // tensor label: I1356
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x3, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, a1, c2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());
    // tensor label: I1357
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x2)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), x2.size());
    dgemm_("T", "N", x3.size()*a1.size()*c2.size(), x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x3.size()*a1.size()*c2.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size(), c2.size(), x2.size());
  out()->put_block(odata, x2, x3, a1, c2);
}

void Task1369::Task_local::compute() {
  const Index a3 = b(0);
  const Index x2 = b(1);
  // tensor label: I1357
  std::unique_ptr<double[]> odata = out()->move_block(a3, x2);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x2);
    sort_indices<0,1,1,1,1,4>(i0data, odata, a3.size(), x2.size());
  }
  out()->put_block(odata, a3, x2);
}

void Task1370::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1320
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), x2.size());
      // tensor label: I1483
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, x0, x1)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x0.size(), x1.size());
      dgemm_("T", "N", x3.size()*x2.size(), x0.size()*x1.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<3,2,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x0.size(), x1.size());
  out()->put_block(odata, x1, x0, x3, x2);
}

void Task1371::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1483
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, x0, x1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
    dscal_(x1.size()*c2.size()*a1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<2,1,0,3,1,1,-1,4>(i0data, odata, x0.size(), a1.size(), c2.size(), x1.size());
  }
  out()->put_block(odata, c2, a1, x0, x1);
}

void Task1372::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1483
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x0, x1), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size());
    // tensor label: I1484
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

void Task1373::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1484
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,1,4>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task1374::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1320
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c3, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c3, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c3.size(), x1.size());
      // tensor label: I1683
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a1, x3, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a1, x3, c3)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), a1.size(), x3.size(), c3.size());
      dgemm_("T", "N", x0.size()*x1.size(), x2.size()*x3.size(), a1.size()*c3.size(),
             1.0, i0data_sorted, a1.size()*c3.size(), i1data_sorted, a1.size()*c3.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size());
  out()->put_block(odata, x1, x0, x3, x2);
}

void Task1375::Task_local::compute() {
  const Index x2 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index c3 = b(3);
  // tensor label: I1683
  std::unique_ptr<double[]> odata = out()->move_block(x2, a1, x3, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a1, x3, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a1, x3, c3), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, c3)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), c3.size());
    // tensor label: I1684
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c2, a1, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c2, a1, x3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), c2.size(), a1.size(), x3.size());
    dgemm_("T", "N", c3.size(), x2.size()*a1.size()*x3.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, c3.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c3.size(), x2.size(), a1.size(), x3.size());
  out()->put_block(odata, x2, a1, x3, c3);
}

void Task1376::Task_local::compute() {
  const Index x2 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x3 = b(3);
  // tensor label: I1684
  std::unique_ptr<double[]> odata = out()->move_block(x2, c2, a1, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, x2);
    sort_indices<3,2,1,0,1,1,-1,4>(i0data, odata, x3.size(), a1.size(), c2.size(), x2.size());
  }
  out()->put_block(odata, x2, c2, a1, x3);
}

void Task1377::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1320
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, c2, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c2.size(), x1.size());
      // tensor label: I1687
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c2, x3, a3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c2, x3, a3)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), c2.size(), x3.size(), a3.size());
      dgemm_("T", "N", x0.size()*x1.size(), x2.size()*x3.size(), c2.size()*a3.size(),
             1.0, i0data_sorted, c2.size()*a3.size(), i1data_sorted, c2.size()*a3.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size());
  out()->put_block(odata, x1, x0, x3, x2);
}

void Task1378::Task_local::compute() {
  const Index x2 = b(0);
  const Index c2 = b(1);
  const Index x3 = b(2);
  const Index a3 = b(3);
  // tensor label: I1687
  std::unique_ptr<double[]> odata = out()->move_block(x2, c2, x3, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, c2, x3, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, c2, x3, a3), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, a1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a3.size(), a1.size());
    // tensor label: I1688
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c2, a1, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c2, a1, x3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), c2.size(), a1.size(), x3.size());
    dgemm_("T", "N", a3.size(), x2.size()*c2.size()*x3.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, a3.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a3.size(), x2.size(), c2.size(), x3.size());
  out()->put_block(odata, x2, c2, x3, a3);
}

void Task1379::Task_local::compute() {
  const Index x2 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x3 = b(3);
  // tensor label: I1688
  std::unique_ptr<double[]> odata = out()->move_block(x2, c2, a1, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, x2);
    sort_indices<3,2,1,0,1,1,1,4>(i0data, odata, x3.size(), a1.size(), c2.size(), x2.size());
  }
  out()->put_block(odata, x2, c2, a1, x3);
}

void Task1380::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1320
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, x2)]);
      sort_indices<2,1,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), x2.size());
      // tensor label: I1718
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, a1, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, a1, c2)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), a1.size(), c2.size());
      dgemm_("T", "N", x2.size()*x3.size(), x1.size()*x0.size(), a1.size()*c2.size(),
             1.0, i0data_sorted, a1.size()*c2.size(), i1data_sorted, a1.size()*c2.size(),
             1.0, odata_sorted, x2.size()*x3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, x3, x2);
}

void Task1381::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a1 = b(2);
  const Index c2 = b(3);
  // tensor label: I1718
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a1, c2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
    // tensor label: I1719
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x1)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), x1.size());
    dgemm_("T", "N", x0.size()*a1.size()*c2.size(), x1.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x0.size()*a1.size()*c2.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), a1.size(), c2.size(), x1.size());
  out()->put_block(odata, x1, x0, a1, c2);
}

void Task1382::Task_local::compute() {
  const Index a3 = b(0);
  const Index x1 = b(1);
  // tensor label: I1719
  std::unique_ptr<double[]> odata = out()->move_block(a3, x1);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x1);
    sort_indices<0,1,1,1,1,4>(i0data, odata, a3.size(), x1.size());
  }
  out()->put_block(odata, a3, x1);
}

void Task1383::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1320
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), x1.size());
      // tensor label: I1845
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, x3, x2)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x3.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), x3.size()*x2.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x3.size(), x2.size());
  out()->put_block(odata, x1, x0, x3, x2);
}

void Task1384::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1845
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, x3, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, x2);
    dscal_(x2.size()*c2.size()*a1.size()*x3.size(), e0_, i0data.get(), 1);
    sort_indices<2,1,0,3,1,1,-1,4>(i0data, odata, x3.size(), a1.size(), c2.size(), x2.size());
  }
  out()->put_block(odata, c2, a1, x3, x2);
}

void Task1385::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1845
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x3, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x2.size(), a3.size());
    // tensor label: I1846
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x3.size());
    dgemm_("T", "N", x2.size(), c2.size()*a1.size()*x3.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x2.size(), c2.size(), a1.size(), x3.size());
  out()->put_block(odata, c2, a1, x3, x2);
}

void Task1386::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x3 = b(3);
  // tensor label: I1846
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,1,4>(i0data, odata, x3.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x3);
}

void Task1387::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1320
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), x2.size());
      // tensor label: I2015
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c2, a1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c2, a1, x0)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), c2.size(), a1.size(), x0.size());
      dgemm_("T", "N", x3.size()*x2.size(), x1.size()*x0.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, x3, x2);
}

void Task1388::Task_local::compute() {
  const Index x1 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I2015
  std::unique_ptr<double[]> odata = out()->move_block(x1, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, x0.size(), a1.size(), c2.size(), x1.size());
  }
  out()->put_block(odata, x1, c2, a1, x0);
}

void Task1389::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1320
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), x1.size());
      // tensor label: I2069
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c2, a1, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c2, a1, x3)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), c2.size(), a1.size(), x3.size());
      dgemm_("T", "N", x0.size()*x1.size(), x2.size()*x3.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size());
  out()->put_block(odata, x1, x0, x3, x2);
}

void Task1390::Task_local::compute() {
  const Index x2 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x3 = b(3);
  // tensor label: I2069
  std::unique_ptr<double[]> odata = out()->move_block(x2, c2, a1, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, x2);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, x3.size(), a1.size(), c2.size(), x2.size());
  }
  out()->put_block(odata, x2, c2, a1, x3);
}

void Task1391::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        for (auto& x0 : *range_[1]) {
          // tensor label: Gamma412
          std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x4, x1, x0);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x5, x4, x1, x0)]);
          sort_indices<1,2,3,4,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x5.size(), x4.size(), x1.size(), x0.size());
          // tensor label: I1328
          std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, x5, x4);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, x5, x4)]);
          sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), x5.size(), x4.size());
          dgemm_("T", "N", ci0.size(), 1, x1.size()*x0.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x1.size()*x0.size()*x5.size()*x4.size(), i1data_sorted, x1.size()*x0.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, ci0.size());
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1392::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  // tensor label: I1328
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x5, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x5, x4), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x5, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x5, x4)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x5.size(), x4.size());
      // tensor label: I1329
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c2, a1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c2, a1, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), c2.size(), a1.size(), x0.size());
      dgemm_("T", "N", x5.size()*x4.size(), x1.size()*x0.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, x5.size()*x4.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, x5, x4);
}

void Task1393::Task_local::compute() {
  const Index x1 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1329
  std::unique_ptr<double[]> odata = out()->move_block(x1, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
    sort_indices<3,2,1,0,1,1,-1,4>(i0data, odata, x0.size(), a1.size(), c2.size(), x1.size());
  }
  out()->put_block(odata, x1, c2, a1, x0);
}

void Task1394::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  // tensor label: I1328
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x5, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x5, x4), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, c1, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, c1, x4)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), c1.size(), x4.size());
      // tensor label: I1372
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, a2, c1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), a2.size(), c1.size());
      dgemm_("T", "N", x5.size()*x4.size(), x1.size()*x0.size(), a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, x5.size()*x4.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, x5, x4);
}

void Task1395::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1372
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
    sort_indices<3,2,1,0,1,1,-1,4>(i0data, odata, c1.size(), a2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x1, x0, a2, c1);
}

void Task1396::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  // tensor label: I1328
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x5, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x5, x4), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x5, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x5, x4)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x5.size(), x4.size());
      // tensor label: I1383
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, a2, c1)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), a2.size(), c1.size());
      dgemm_("T", "N", x5.size()*x4.size(), x1.size()*x0.size(), a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, x5.size()*x4.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, x5, x4);
}

void Task1397::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1383
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, c1.size(), a2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x1, x0, a2, c1);
}

void Task1398::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  // tensor label: I1328
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x5, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x5, x4), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x0.size(), x1.size());
      // tensor label: I1691
      std::unique_ptr<double[]> i1data = in(1)->get_block(x4, c2, a1, x5);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, c2, a1, x5)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x4.size(), c2.size(), a1.size(), x5.size());
      dgemm_("T", "N", x0.size()*x1.size(), x4.size()*x5.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x4.size(), x5.size());
  out()->put_block(odata, x1, x0, x5, x4);
}

void Task1399::Task_local::compute() {
  const Index x4 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x5 = b(3);
  // tensor label: I1691
  std::unique_ptr<double[]> odata = out()->move_block(x4, c2, a1, x5);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c2, x4);
    sort_indices<3,2,1,0,1,1,-1,4>(i0data, odata, x5.size(), a1.size(), c2.size(), x4.size());
  }
  out()->put_block(odata, x4, c2, a1, x5);
}

