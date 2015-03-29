//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks24.cc
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

#include <src/smith/MRCI_tasks24.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task1150::Task_local::compute() {
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
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, c4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, c4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), c4.size());
      // tensor label: I1392
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c4, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c4, x1, x0)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c4.size(), x1.size(), x0.size());
      dgemm_("T", "N", a1.size()*c2.size(), a3.size()*x0.size(), c4.size()*x1.size(),
             1.0, i0data_sorted, c4.size()*x1.size(), i1data_sorted, c4.size()*x1.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<2,1,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1151::Task_local::compute() {
  const Index a3 = b(0);
  const Index c4 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1392
  std::unique_ptr<double[]> odata = out()->move_block(a3, c4, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c4, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c4, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma29
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I1393
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a3, c4, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a3, c4, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a3.size(), c4.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), a3.size()*c4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a3.size(), c4.size());
  out()->put_block(odata, a3, c4, x1, x0);
}

void Task1152::Task_local::compute() {
  const Index x3 = b(0);
  const Index a3 = b(1);
  const Index c4 = b(2);
  const Index x2 = b(3);
  // tensor label: I1393
  std::unique_ptr<double[]> odata = out()->move_block(x3, a3, c4, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c4, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), a3.size(), c4.size(), x2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c4, a3, x3, x2);
    sort_indices<2,1,0,3,1,1,-2,1>(i1data, odata, c4.size(), a3.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, x3, a3, c4, x2);
}

void Task1153::Task_local::compute() {
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
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, c4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, c4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), c4.size());
      // tensor label: I1395
      std::unique_ptr<double[]> i1data = in(1)->get_block(a1, c4, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, c4, x0, x1)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, a1.size(), c4.size(), x0.size(), x1.size());
      dgemm_("T", "N", a3.size()*c2.size(), a1.size()*x0.size(), c4.size()*x1.size(),
             1.0, i0data_sorted, c4.size()*x1.size(), i1data_sorted, c4.size()*x1.size(),
             1.0, odata_sorted, a3.size()*c2.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1154::Task_local::compute() {
  const Index a1 = b(0);
  const Index c4 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1395
  std::unique_ptr<double[]> odata = out()->move_block(a1, c4, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, c4, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, c4, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma27
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      // tensor label: I1396
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

void Task1155::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index c4 = b(2);
  const Index x2 = b(3);
  // tensor label: I1396
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1, c4, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c4, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), a1.size(), c4.size(), x2.size());
  }
  out()->put_block(odata, x3, a1, c4, x2);
}

void Task1156::Task_local::compute() {
  const Index a1 = b(0);
  const Index c4 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1395
  std::unique_ptr<double[]> odata = out()->move_block(a1, c4, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, c4, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, c4, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma29
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I1426
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

void Task1157::Task_local::compute() {
  const Index c4 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1426
  std::unique_ptr<double[]> odata = out()->move_block(c4, a1, x3, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a1, x3, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c4.size(), a1.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, c4, a1, x3, x2);
}

void Task1158::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, a4, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, a4, a1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), a4.size(), a1.size());
      // tensor label: I1398
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c2, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c2, x0, x1)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a4.size(), c2.size(), x0.size(), x1.size());
      dgemm_("T", "N", a3.size()*a1.size(), c2.size()*x0.size(), a4.size()*x1.size(),
             1.0, i0data_sorted, a4.size()*x1.size(), i1data_sorted, a4.size()*x1.size(),
             1.0, odata_sorted, a3.size()*a1.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, a3.size(), a1.size(), c2.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1159::Task_local::compute() {
  const Index a4 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1398
  std::unique_ptr<double[]> odata = out()->move_block(a4, c2, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma27
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      // tensor label: I1399
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a4, c2, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a4, c2, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a4.size(), c2.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), a4.size()*c2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a4.size(), c2.size());
  out()->put_block(odata, a4, c2, x0, x1);
}

void Task1160::Task_local::compute() {
  const Index x3 = b(0);
  const Index a4 = b(1);
  const Index c2 = b(2);
  const Index x2 = b(3);
  // tensor label: I1399
  std::unique_ptr<double[]> odata = out()->move_block(x3, a4, c2, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c2, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), a4.size(), c2.size(), x2.size());
  }
  out()->put_block(odata, x3, a4, c2, x2);
}

void Task1161::Task_local::compute() {
  const Index a4 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1398
  std::unique_ptr<double[]> odata = out()->move_block(a4, c2, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma29
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I1429
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a4, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a4, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), a4.size(), x3.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), c2.size()*a4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a4.size());
  out()->put_block(odata, a4, c2, x0, x1);
}

void Task1162::Task_local::compute() {
  const Index c2 = b(0);
  const Index a4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1429
  std::unique_ptr<double[]> odata = out()->move_block(c2, a4, x3, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a4, x3, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c2.size(), a4.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, c2, a4, x3, x2);
}

void Task1163::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, a4, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, a4, a3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), a4.size(), a3.size());
      // tensor label: I1401
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c2, x1, x0)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), c2.size(), x1.size(), x0.size());
      dgemm_("T", "N", a1.size()*a3.size(), c2.size()*x0.size(), a4.size()*x1.size(),
             1.0, i0data_sorted, a4.size()*x1.size(), i1data_sorted, a4.size()*x1.size(),
             1.0, odata_sorted, a1.size()*a3.size());
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), c2.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1164::Task_local::compute() {
  const Index a4 = b(0);
  const Index c2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1401
  std::unique_ptr<double[]> odata = out()->move_block(a4, c2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c2, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma29
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I1402
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a4, c2, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a4, c2, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a4.size(), c2.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), a4.size()*c2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a4.size(), c2.size());
  out()->put_block(odata, a4, c2, x1, x0);
}

void Task1165::Task_local::compute() {
  const Index x3 = b(0);
  const Index a4 = b(1);
  const Index c2 = b(2);
  const Index x2 = b(3);
  // tensor label: I1402
  std::unique_ptr<double[]> odata = out()->move_block(x3, a4, c2, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c2, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), a4.size(), c2.size(), x2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c2, a4, x3, x2);
    sort_indices<2,1,0,3,1,1,2,1>(i1data, odata, c2.size(), a4.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, x3, a4, c2, x2);
}

void Task1166::Task_local::compute() {
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
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, x5, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, x5, x4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), x5.size(), x4.size());
      // tensor label: I1404
      std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x5, x4, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x5, x4, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, a1.size(), x5.size(), x4.size(), x0.size());
      dgemm_("T", "N", c2.size()*a3.size(), a1.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), a1.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1167::Task_local::compute() {
  const Index a1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  // tensor label: I1404
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
        // tensor label: I1405
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

void Task1168::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1405
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1, x2, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), a1.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, a1, x2, x1);
}

void Task1169::Task_local::compute() {
  const Index a1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  // tensor label: I1404
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
        // tensor label: I1411
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

void Task1170::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  // tensor label: I1411
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, a1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, a1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x2.size(), x1.size(), a1.size());
  }
  out()->put_block(odata, x3, x2, x1, a1);
}

void Task1171::Task_local::compute() {
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
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x5, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x5, x4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x5.size(), x4.size());
      // tensor label: I1407
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x5, x4, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x5, x4, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), x5.size(), x4.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), a3.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1172::Task_local::compute() {
  const Index a3 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  // tensor label: I1407
  std::unique_ptr<double[]> odata = out()->move_block(a3, x5, x4, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x5, x4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x5, x4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma49
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x0, x2, x1)]);
        sort_indices<2,4,5,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: I1408
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a3, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), a3.size(), x2.size(), x1.size());
        dgemm_("T", "N", x5.size()*x4.size()*x0.size(), a3.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x4.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x0.size(), a3.size());
  out()->put_block(odata, a3, x5, x4, x0);
}

void Task1173::Task_local::compute() {
  const Index x3 = b(0);
  const Index a3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1408
  std::unique_ptr<double[]> odata = out()->move_block(x3, a3, x2, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x3.size(), a3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, a3, x2, x1);
}

void Task1174::Task_local::compute() {
  const Index a3 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  // tensor label: I1407
  std::unique_ptr<double[]> odata = out()->move_block(a3, x5, x4, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x5, x4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x5, x4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma240
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x2, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x2, x1, x0)]);
        sort_indices<2,3,4,0,1,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());
        // tensor label: I1414
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, a3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, a3)]);
        sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), a3.size());
        dgemm_("T", "N", x5.size()*x4.size()*x0.size(), a3.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x4.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x0.size(), a3.size());
  out()->put_block(odata, a3, x5, x4, x0);
}

void Task1175::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a3 = b(3);
  // tensor label: I1414
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, a3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, a3);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x3.size(), x2.size(), x1.size(), a3.size());
  }
  out()->put_block(odata, x3, x2, x1, a3);
}

void Task1176::Task_local::compute() {
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
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x2.size(), x1.size());
      // tensor label: I1434
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x0, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), x0.size(), x2.size(), x1.size());
      dgemm_("T", "N", c2.size()*a1.size(), a3.size()*x0.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1177::Task_local::compute() {
  const Index a3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1434
  std::unique_ptr<double[]> odata = out()->move_block(a3, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma50
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: I1435
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a3, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a3, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x5.size(), a3.size(), x4.size(), x3.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a3.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a3.size());
  out()->put_block(odata, a3, x0, x2, x1);
}

void Task1178::Task_local::compute() {
  const Index x5 = b(0);
  const Index a3 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I1435
  std::unique_ptr<double[]> odata = out()->move_block(x5, a3, x4, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, x4, x3);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x5.size(), a3.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, x5, a3, x4, x3);
}

void Task1179::Task_local::compute() {
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
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), x2.size(), x1.size());
      // tensor label: I1437
      std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x0, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a1.size(), x0.size(), x2.size(), x1.size());
      dgemm_("T", "N", c2.size()*a3.size(), a1.size()*x0.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), a1.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1180::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1437
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma50
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: I1438
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a1, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a1, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->put_block(odata, a1, x0, x2, x1);
}

void Task1181::Task_local::compute() {
  const Index x5 = b(0);
  const Index a1 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I1438
  std::unique_ptr<double[]> odata = out()->move_block(x5, a1, x4, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x5.size(), a1.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, x5, a1, x4, x3);
}

void Task1182::Task_local::compute() {
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
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x2, x1, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x2, x1, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), x2.size(), x1.size(), a1.size());
      // tensor label: I1440
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x2, x1, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), a3.size()*x0.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1183::Task_local::compute() {
  const Index a3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1440
  std::unique_ptr<double[]> odata = out()->move_block(a3, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma252
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2, x4, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x2, x4, x3, x1, x0)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
        // tensor label: I1441
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a3, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a3, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x5.size(), a3.size(), x4.size(), x3.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), a3.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), a3.size());
  out()->put_block(odata, a3, x2, x1, x0);
}

void Task1184::Task_local::compute() {
  const Index x5 = b(0);
  const Index a3 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I1441
  std::unique_ptr<double[]> odata = out()->move_block(x5, a3, x4, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, x4, x3);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x5.size(), a3.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, x5, a3, x4, x3);
}

void Task1185::Task_local::compute() {
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
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x2, x1, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x2, x1, a3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), x2.size(), x1.size(), a3.size());
      // tensor label: I1443
      std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x0, x1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x0, x1, x2)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, a1.size(), x0.size(), x1.size(), x2.size());
      dgemm_("T", "N", c2.size()*a3.size(), a1.size()*x0.size(), x1.size()*x2.size(),
             1.0, i0data_sorted, x1.size()*x2.size(), i1data_sorted, x1.size()*x2.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), a1.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1186::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: I1443
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x1, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, x2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma31
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x1, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x1, x2)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x1.size(), x2.size());
        // tensor label: I1444
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a1, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a1, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
        dgemm_("T", "N", x0.size()*x1.size()*x2.size(), a1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x1.size()*x2.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), a1.size());
  out()->put_block(odata, a1, x0, x1, x2);
}

void Task1187::Task_local::compute() {
  const Index x5 = b(0);
  const Index a1 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I1444
  std::unique_ptr<double[]> odata = out()->move_block(x5, a1, x4, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x5.size(), a1.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, x5, a1, x4, x3);
}

void Task1188::Task_local::compute() {
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
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a1, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a1, c2, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x2.size(), a1.size(), c2.size(), x1.size());
      // tensor label: I1446
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x1, x2, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x1, x2, x0)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), x1.size(), x2.size(), x0.size());
      dgemm_("T", "N", a1.size()*c2.size(), a3.size()*x0.size(), x1.size()*x2.size(),
             1.0, i0data_sorted, x1.size()*x2.size(), i1data_sorted, x1.size()*x2.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<2,1,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1189::Task_local::compute() {
  const Index a3 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index x0 = b(3);
  // tensor label: I1446
  std::unique_ptr<double[]> odata = out()->move_block(a3, x1, x2, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x1, x2, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x1, x2, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma471
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x1, x4, x3, x2, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x1, x4, x3, x2, x0)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x1.size(), x4.size(), x3.size(), x2.size(), x0.size());
        // tensor label: I1447
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a3, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a3, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x5.size(), a3.size(), x4.size(), x3.size());
        dgemm_("T", "N", x1.size()*x2.size()*x0.size(), a3.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x1.size()*x2.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x2.size(), x0.size(), a3.size());
  out()->put_block(odata, a3, x1, x2, x0);
}

void Task1190::Task_local::compute() {
  const Index x5 = b(0);
  const Index a3 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I1447
  std::unique_ptr<double[]> odata = out()->move_block(x5, a3, x4, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, x4, x3);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x5.size(), a3.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, x5, a3, x4, x3);
}

void Task1191::Task_local::compute() {
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
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a3, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a3, c2, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x2.size(), a3.size(), c2.size(), x1.size());
      // tensor label: I1449
      std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x0, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a1.size(), x0.size(), x2.size(), x1.size());
      dgemm_("T", "N", a3.size()*c2.size(), a1.size()*x0.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, a3.size()*c2.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1192::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1449
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma50
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: I1450
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a1, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a1, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->put_block(odata, a1, x0, x2, x1);
}

void Task1193::Task_local::compute() {
  const Index x5 = b(0);
  const Index a1 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I1450
  std::unique_ptr<double[]> odata = out()->move_block(x5, a1, x4, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x5.size(), a1.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, x5, a1, x4, x3);
}

void Task1194::Task_local::compute() {
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
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x1, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x1, c2, a1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x2.size(), x1.size(), c2.size(), a1.size());
      // tensor label: I1452
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x0, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), x0.size(), x2.size(), x1.size());
      dgemm_("T", "N", c2.size()*a1.size(), a3.size()*x0.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1195::Task_local::compute() {
  const Index a3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1452
  std::unique_ptr<double[]> odata = out()->move_block(a3, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma50
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: I1453
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a3, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a3, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x5.size(), a3.size(), x4.size(), x3.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a3.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a3.size());
  out()->put_block(odata, a3, x0, x2, x1);
}

void Task1196::Task_local::compute() {
  const Index x5 = b(0);
  const Index a3 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I1453
  std::unique_ptr<double[]> odata = out()->move_block(x5, a3, x4, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, x4, x3);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x5.size(), a3.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, x5, a3, x4, x3);
}

void Task1197::Task_local::compute() {
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
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x1, c2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x1, c2, a3)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x2.size(), x1.size(), c2.size(), a3.size());
      // tensor label: I1455
      std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x0, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a1.size(), x0.size(), x2.size(), x1.size());
      dgemm_("T", "N", c2.size()*a3.size(), a1.size()*x0.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), a1.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1198::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1455
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma50
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: I1456
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a1, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a1, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->put_block(odata, a1, x0, x2, x1);
}

void Task1199::Task_local::compute() {
  const Index x5 = b(0);
  const Index a1 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I1456
  std::unique_ptr<double[]> odata = out()->move_block(x5, a1, x4, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x5.size(), a1.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, x5, a1, x4, x3);
}

#endif
