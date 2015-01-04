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
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1255
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3, a2, c1);
}

void Task1251::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  // tensor label: I1251
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c1 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c1, x0)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c1.size(), x0.size());
        // tensor label: I1614
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c3, a2, c1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c3, a2, c1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, x3.size(), c3.size(), a2.size(), c1.size());
        dgemm_("T", "N", x0.size(), x3.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x3.size());
  out()->put_block(odata, x0, x3);
}

void Task1252::Task_local::compute() {
  const Index x3 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1614
  std::unique_ptr<double[]> odata = out()->move_block(x3, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x3);
    sort_indices<3,2,1,0,1,1,-1,4>(i0data, odata, c1.size(), a2.size(), c3.size(), x3.size());
  }
  out()->put_block(odata, x3, c3, a2, c1);
}

void Task1253::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  // tensor label: I1251
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        // tensor label: I1617
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c3, a2, c1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c3, a2, c1)]);
        sort_indices<3,2,1,0,0,1,1,1>(i1data, i1data_sorted, x3.size(), c3.size(), a2.size(), c1.size());
        dgemm_("T", "N", x0.size(), x3.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x3.size());
  out()->put_block(odata, x0, x3);
}

void Task1254::Task_local::compute() {
  const Index x3 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1617
  std::unique_ptr<double[]> odata = out()->move_block(x3, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x3);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, c1.size(), a2.size(), c3.size(), x3.size());
  }
  out()->put_block(odata, x3, c3, a2, c1);
}

void Task1255::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma394
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x0, x1)]);
      sort_indices<1,2,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x0.size(), x1.size());
      // tensor label: I1257
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
      dgemm_("T", "N", ci0.size(), 1, x0.size()*x1.size(),
             1.0, i0data_sorted, x0.size()*x1.size(), i1data_sorted, x0.size()*x1.size(),
             1.0, odata_sorted, ci0.size());
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1256::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I1257
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c3, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a2, c3, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c3.size(), x1.size());
        // tensor label: I1258
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c4)]);
        sort_indices<3,2,1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c4.size());
        dgemm_("T", "N", x1.size(), x0.size(), c3.size()*a2.size()*c4.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c4.size(), i1data_sorted, c3.size()*a2.size()*c4.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x0, x1);
}

void Task1257::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c4 = b(3);
  // tensor label: I1258
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3, a2, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3, a2, c4), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, c4)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), c4.size());
    // tensor label: I1259
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c1.size());
    dgemm_("T", "N", c4.size(), x0.size()*c3.size()*a2.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, c4.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c4.size(), x0.size(), c3.size(), a2.size());
  out()->put_block(odata, x0, c3, a2, c4);
}

void Task1258::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1259
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3, a2, c1);
}

void Task1259::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I1257
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c4 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c4, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c4, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c4.size(), x1.size());
        // tensor label: I1262
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c4)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c4.size());
        dgemm_("T", "N", x1.size(), x0.size(), c3.size()*a2.size()*c4.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c4.size(), i1data_sorted, c3.size()*a2.size()*c4.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x0, x1);
}

void Task1260::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c4 = b(3);
  // tensor label: I1262
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3, a2, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3, a2, c4), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, c4)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), c4.size());
    // tensor label: I1263
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c1.size());
    dgemm_("T", "N", c4.size(), x0.size()*c3.size()*a2.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, c4.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c4.size(), x0.size(), c3.size(), a2.size());
  out()->put_block(odata, x0, c3, a2, c4);
}

void Task1261::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1263
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<3,2,1,0,1,1,1,4>(i0data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3, a2, c1);
}

void Task1262::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I1257
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c1 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c1, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a2, c1, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c1.size(), x1.size());
        // tensor label: I1266
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a2, c1, c4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a2, c1, c4)]);
        sort_indices<3,1,2,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), a2.size(), c1.size(), c4.size());
        dgemm_("T", "N", x1.size(), x0.size(), a2.size()*c1.size()*c4.size(),
               1.0, i0data_sorted, a2.size()*c1.size()*c4.size(), i1data_sorted, a2.size()*c1.size()*c4.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x0, x1);
}

void Task1263::Task_local::compute() {
  const Index x0 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index c4 = b(3);
  // tensor label: I1266
  std::unique_ptr<double[]> odata = out()->move_block(x0, a2, c1, c4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a2, c1, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a2, c1, c4), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, c4)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), c4.size());
    // tensor label: I1267
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c1.size());
    dgemm_("T", "N", c4.size(), x0.size()*a2.size()*c1.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, c4.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c4.size(), x0.size(), a2.size(), c1.size());
  out()->put_block(odata, x0, a2, c1, c4);
}

void Task1264::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1267
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<3,2,1,0,1,1,1,4>(i0data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3, a2, c1);
}

void Task1265::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I1257
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      for (auto& c1 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, c1, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4, c1, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), c1.size(), x1.size());
        // tensor label: I1270
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, c1, a4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, c1, a4)]);
        sort_indices<1,3,2,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), c1.size(), a4.size());
        dgemm_("T", "N", x1.size(), x0.size(), c3.size()*c1.size()*a4.size(),
               1.0, i0data_sorted, c3.size()*c1.size()*a4.size(), i1data_sorted, c3.size()*c1.size()*a4.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x0, x1);
}

void Task1266::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index c1 = b(2);
  const Index a4 = b(3);
  // tensor label: I1270
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, c1, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3, c1, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3, c1, a4), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a4.size(), a2.size());
    // tensor label: I1271
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c1.size());
    dgemm_("T", "N", a4.size(), x0.size()*c3.size()*c1.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, a4.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a4.size(), x0.size(), c3.size(), c1.size());
  out()->put_block(odata, x0, c3, c1, a4);
}

void Task1267::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1271
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<3,2,1,0,1,1,-1,4>(i0data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3, a2, c1);
}

void Task1268::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I1257
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c4 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c4, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c4, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c4.size(), x1.size());
        // tensor label: I1274
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a2, c1, c4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a2, c1, c4)]);
        sort_indices<2,1,3,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), a2.size(), c1.size(), c4.size());
        dgemm_("T", "N", x1.size(), x0.size(), a2.size()*c1.size()*c4.size(),
               1.0, i0data_sorted, a2.size()*c1.size()*c4.size(), i1data_sorted, a2.size()*c1.size()*c4.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x0, x1);
}

void Task1269::Task_local::compute() {
  const Index x0 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index c4 = b(3);
  // tensor label: I1274
  std::unique_ptr<double[]> odata = out()->move_block(x0, a2, c1, c4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a2, c1, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a2, c1, c4), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, c4)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), c4.size());
    // tensor label: I1275
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c1.size());
    dgemm_("T", "N", c4.size(), x0.size()*a2.size()*c1.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, c4.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c4.size(), x0.size(), a2.size(), c1.size());
  out()->put_block(odata, x0, a2, c1, c4);
}

void Task1270::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1275
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3, a2, c1);
}

void Task1271::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I1257
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x1.size());
        // tensor label: I1278
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, c1, a4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, c1, a4)]);
        sort_indices<2,3,1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), c1.size(), a4.size());
        dgemm_("T", "N", x1.size(), x0.size(), c3.size()*c1.size()*a4.size(),
               1.0, i0data_sorted, c3.size()*c1.size()*a4.size(), i1data_sorted, c3.size()*c1.size()*a4.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x0, x1);
}

void Task1272::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index c1 = b(2);
  const Index a4 = b(3);
  // tensor label: I1278
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, c1, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3, c1, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3, c1, a4), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a4.size(), a2.size());
    // tensor label: I1279
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c1.size());
    dgemm_("T", "N", a4.size(), x0.size()*c3.size()*c1.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, a4.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a4.size(), x0.size(), c3.size(), c1.size());
  out()->put_block(odata, x0, c3, c1, a4);
}

void Task1273::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1279
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3, a2, c1);
}

void Task1274::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I1257
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), x1.size());
    // tensor label: I1298
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size());
    dgemm_("T", "N", x1.size(), x0.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x0, x1);
}

void Task1275::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  // tensor label: I1298
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& c3 : *range_[0]) {
      for (auto& a2 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
        // tensor label: I1299
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c1)]);
        sort_indices<3,1,2,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c1.size());
        dgemm_("T", "N", a4.size(), x0.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, a4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a4.size(), x0.size());
  out()->put_block(odata, x0, a4);
}

void Task1276::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1299
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3, a2, c1);
}

void Task1277::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  // tensor label: I1298
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
        // tensor label: I1303
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c1)]);
        sort_indices<3,2,1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c1.size());
        dgemm_("T", "N", a4.size(), x0.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, a4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a4.size(), x0.size());
  out()->put_block(odata, x0, a4);
}

void Task1278::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1303
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3, a2, c1);
}

void Task1279::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I1257
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a2.size());
    // tensor label: I1441
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x1)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a2.size(), x1.size());
    dgemm_("T", "N", x0.size(), x1.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->put_block(odata, x0, x1);
}

void Task1280::Task_local::compute() {
  const Index a2 = b(0);
  const Index x1 = b(1);
  // tensor label: I1441
  std::unique_ptr<double[]> odata = out()->move_block(a2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x1.size());
        // tensor label: I1442
        std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3, a2, c1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3, a2, c1)]);
        sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size(), a2.size(), c1.size());
        dgemm_("T", "N", x1.size(), a2.size(), a4.size()*c3.size()*c1.size(),
               1.0, i0data_sorted, a4.size()*c3.size()*c1.size(), i1data_sorted, a4.size()*c3.size()*c1.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), a2.size());
  out()->put_block(odata, a2, x1);
}

void Task1281::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1442
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->put_block(odata, a4, c3, a2, c1);
}

void Task1282::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I1257
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a4)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a4.size());
    // tensor label: I1445
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x1)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), x1.size());
    dgemm_("T", "N", x0.size(), x1.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->put_block(odata, x0, x1);
}

void Task1283::Task_local::compute() {
  const Index a4 = b(0);
  const Index x1 = b(1);
  // tensor label: I1445
  std::unique_ptr<double[]> odata = out()->move_block(a4, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x1.size());
        // tensor label: I1446
        std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3, a2, c1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3, a2, c1)]);
        sort_indices<3,2,1,0,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size(), a2.size(), c1.size());
        dgemm_("T", "N", x1.size(), a4.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), a4.size());
  out()->put_block(odata, a4, x1);
}

void Task1284::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1446
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->put_block(odata, a4, c3, a2, c1);
}

void Task1285::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I1257
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c3, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a2, c3, x0)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c3.size(), x0.size());
        // tensor label: I1620
        std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c3, a2, c4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c3, a2, c4)]);
        sort_indices<3,2,1,0,0,1,1,1>(i1data, i1data_sorted, x1.size(), c3.size(), a2.size(), c4.size());
        dgemm_("T", "N", x0.size(), x1.size(), c3.size()*a2.size()*c4.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c4.size(), i1data_sorted, c3.size()*a2.size()*c4.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->put_block(odata, x0, x1);
}

void Task1286::Task_local::compute() {
  const Index x1 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c4 = b(3);
  // tensor label: I1620
  std::unique_ptr<double[]> odata = out()->move_block(x1, c3, a2, c4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, c3, a2, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, c3, a2, c4), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, c4)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), c4.size());
    // tensor label: I1621
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c3, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c3, a2, c1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), c3.size(), a2.size(), c1.size());
    dgemm_("T", "N", c4.size(), x1.size()*c3.size()*a2.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, c4.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c4.size(), x1.size(), c3.size(), a2.size());
  out()->put_block(odata, x1, c3, a2, c4);
}

void Task1287::Task_local::compute() {
  const Index x1 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1621
  std::unique_ptr<double[]> odata = out()->move_block(x1, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x1);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, c1.size(), a2.size(), c3.size(), x1.size());
  }
  out()->put_block(odata, x1, c3, a2, c1);
}

void Task1288::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I1257
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c4 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c4, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c4, x0)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c4.size(), x0.size());
        // tensor label: I1624
        std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c3, a2, c4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c3, a2, c4)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, x1.size(), c3.size(), a2.size(), c4.size());
        dgemm_("T", "N", x0.size(), x1.size(), c3.size()*a2.size()*c4.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c4.size(), i1data_sorted, c3.size()*a2.size()*c4.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->put_block(odata, x0, x1);
}

void Task1289::Task_local::compute() {
  const Index x1 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c4 = b(3);
  // tensor label: I1624
  std::unique_ptr<double[]> odata = out()->move_block(x1, c3, a2, c4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, c3, a2, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, c3, a2, c4), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, c4)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), c4.size());
    // tensor label: I1625
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c3, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c3, a2, c1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), c3.size(), a2.size(), c1.size());
    dgemm_("T", "N", c4.size(), x1.size()*c3.size()*a2.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, c4.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c4.size(), x1.size(), c3.size(), a2.size());
  out()->put_block(odata, x1, c3, a2, c4);
}

void Task1290::Task_local::compute() {
  const Index x1 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1625
  std::unique_ptr<double[]> odata = out()->move_block(x1, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x1);
    sort_indices<3,2,1,0,1,1,1,4>(i0data, odata, c1.size(), a2.size(), c3.size(), x1.size());
  }
  out()->put_block(odata, x1, c3, a2, c1);
}

void Task1291::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I1257
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c1 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a2, c1, x0)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c1.size(), x0.size());
        // tensor label: I1628
        std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a2, c1, c4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a2, c1, c4)]);
        sort_indices<3,1,2,0,0,1,1,1>(i1data, i1data_sorted, x1.size(), a2.size(), c1.size(), c4.size());
        dgemm_("T", "N", x0.size(), x1.size(), a2.size()*c1.size()*c4.size(),
               1.0, i0data_sorted, a2.size()*c1.size()*c4.size(), i1data_sorted, a2.size()*c1.size()*c4.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->put_block(odata, x0, x1);
}

void Task1292::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index c4 = b(3);
  // tensor label: I1628
  std::unique_ptr<double[]> odata = out()->move_block(x1, a2, c1, c4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, a2, c1, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, a2, c1, c4), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, c4)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), c4.size());
    // tensor label: I1629
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c3, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c3, a2, c1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), c3.size(), a2.size(), c1.size());
    dgemm_("T", "N", c4.size(), x1.size()*a2.size()*c1.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, c4.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c4.size(), x1.size(), a2.size(), c1.size());
  out()->put_block(odata, x1, a2, c1, c4);
}

void Task1293::Task_local::compute() {
  const Index x1 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1629
  std::unique_ptr<double[]> odata = out()->move_block(x1, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x1);
    sort_indices<3,2,1,0,1,1,1,4>(i0data, odata, c1.size(), a2.size(), c3.size(), x1.size());
  }
  out()->put_block(odata, x1, c3, a2, c1);
}

void Task1294::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I1257
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      for (auto& c1 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, c1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4, c1, x0)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), c1.size(), x0.size());
        // tensor label: I1632
        std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c3, c1, a4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c3, c1, a4)]);
        sort_indices<1,3,2,0,0,1,1,1>(i1data, i1data_sorted, x1.size(), c3.size(), c1.size(), a4.size());
        dgemm_("T", "N", x0.size(), x1.size(), c3.size()*c1.size()*a4.size(),
               1.0, i0data_sorted, c3.size()*c1.size()*a4.size(), i1data_sorted, c3.size()*c1.size()*a4.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->put_block(odata, x0, x1);
}

void Task1295::Task_local::compute() {
  const Index x1 = b(0);
  const Index c3 = b(1);
  const Index c1 = b(2);
  const Index a4 = b(3);
  // tensor label: I1632
  std::unique_ptr<double[]> odata = out()->move_block(x1, c3, c1, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, c3, c1, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, c3, c1, a4), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a4.size(), a2.size());
    // tensor label: I1633
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c3, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c3, a2, c1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), c3.size(), a2.size(), c1.size());
    dgemm_("T", "N", a4.size(), x1.size()*c3.size()*c1.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, a4.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a4.size(), x1.size(), c3.size(), c1.size());
  out()->put_block(odata, x1, c3, c1, a4);
}

void Task1296::Task_local::compute() {
  const Index x1 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1633
  std::unique_ptr<double[]> odata = out()->move_block(x1, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x1);
    sort_indices<3,2,1,0,1,1,-1,4>(i0data, odata, c1.size(), a2.size(), c3.size(), x1.size());
  }
  out()->put_block(odata, x1, c3, a2, c1);
}

void Task1297::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I1257
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c4 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c4, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c4, x0)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c4.size(), x0.size());
        // tensor label: I1636
        std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a2, c1, c4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a2, c1, c4)]);
        sort_indices<2,1,3,0,0,1,1,1>(i1data, i1data_sorted, x1.size(), a2.size(), c1.size(), c4.size());
        dgemm_("T", "N", x0.size(), x1.size(), a2.size()*c1.size()*c4.size(),
               1.0, i0data_sorted, a2.size()*c1.size()*c4.size(), i1data_sorted, a2.size()*c1.size()*c4.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->put_block(odata, x0, x1);
}

void Task1298::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index c4 = b(3);
  // tensor label: I1636
  std::unique_ptr<double[]> odata = out()->move_block(x1, a2, c1, c4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, a2, c1, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, a2, c1, c4), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, c4)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), c4.size());
    // tensor label: I1637
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c3, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c3, a2, c1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), c3.size(), a2.size(), c1.size());
    dgemm_("T", "N", c4.size(), x1.size()*a2.size()*c1.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, c4.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c4.size(), x1.size(), a2.size(), c1.size());
  out()->put_block(odata, x1, a2, c1, c4);
}

void Task1299::Task_local::compute() {
  const Index x1 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1637
  std::unique_ptr<double[]> odata = out()->move_block(x1, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x1);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, c1.size(), a2.size(), c3.size(), x1.size());
  }
  out()->put_block(odata, x1, c3, a2, c1);
}

