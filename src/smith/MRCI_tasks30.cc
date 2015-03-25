//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks30.cc
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

#include <src/smith/MRCI_tasks30.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task1450::Task_local::compute() {
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I1798
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma33
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I1801
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a1, c2, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), a1.size(), c2.size(), a3.size());
    dgemm_("T", "N", x0.size(), a1.size()*c2.size()*a3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, x0.size(), a1.size(), c2.size(), a3.size());
  out()->put_block(odata, x0, a3, c2, a1);
}

void Task1451::Task_local::compute() {
  const Index x1 = b(0);
  const Index a1 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  // tensor label: I1801
  std::unique_ptr<double[]> odata = out()->move_block(x1, a1, c2, a3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x1.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, x1, a1, c2, a3);
}

void Task1452::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(x1, a2, x0, a1);
  {
    // tensor label: I1802
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, a1, a2);
    sort_indices<1,3,0,2,1,1,1,1>(i0data, odata, x0.size(), x1.size(), a1.size(), a2.size());
  }
  out()->put_block(odata, x1, a2, x0, a1);
}

void Task1453::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I1802
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());
      // tensor label: I1803
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      dgemm_("T", "N", a1.size()*a2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, a1, a2);
}

void Task1454::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1803
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma53
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}

#endif
