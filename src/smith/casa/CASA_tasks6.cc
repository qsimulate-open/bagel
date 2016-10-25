//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASA_tasks6.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#include <src/smith/casa/CASA_tasks6.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASA;

void Task250::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I340
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, a1, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c2, a1, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x1, x0), 0.0);
  // tensor label: Gamma29
  std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
  // tensor label: h1
  std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1)]);
  sort_indices<0,1,0,1,-1,1>(i1data, i1data_sorted, c2.size(), a1.size());
  dgemm_("T", "N", x1.size()*x0.size(), c2.size()*a1.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, x1.size()*x0.size());
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a1.size());
  out()->add_block(odata, c2, a1, x1, x0);
}

void Task251::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1, c1, a2), 0.0);
  {
    // tensor label: I348
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x1, x0);
    sort_indices<3,2,0,1,1,1,1,1>(i0data, odata, c1.size(), a2.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x0, x1, c1, a2);
}

void Task252::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I348
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma25
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I349
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), x3.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), c1.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), a2.size());
  out()->add_block(odata, c1, a2, x1, x0);
}

void Task253::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I349
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, x3, x2), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x3, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c1.size(), a2.size(), x3.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x3, a2, c1, x2);
    sort_indices<2,1,0,3,1,1,-1,2>(i1data, odata, x3.size(), a2.size(), c1.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i2data = in(0)->get_block(x3, x2, c1, a2);
    sort_indices<2,3,0,1,1,1,1,1>(i2data, odata, x3.size(), x2.size(), c1.size(), a2.size());
  }
  out()->add_block(odata, c1, a2, x3, x2);
}

void Task254::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I348
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, x1, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: Gamma5
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x3, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3, x2, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x3, x2, a2)]);
      sort_indices<2,1,0,3,0,1,1,2>(i1data, i1data_sorted, c1.size(), x3.size(), x2.size(), a2.size());
      dgemm_("T", "N", x1.size()*x0.size(), c1.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), a2.size());
  out()->add_block(odata, c1, a2, x1, x0);
}

void Task255::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I348
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, x1, x0), 0.0);
  // tensor label: Gamma29
  std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
  // tensor label: h1
  std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2)]);
  sort_indices<0,1,0,1,2,1>(i1data, i1data_sorted, c1.size(), a2.size());
  dgemm_("T", "N", x1.size()*x0.size(), c1.size()*a2.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, x1.size()*x0.size());
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), a2.size());
  out()->add_block(odata, c1, a2, x1, x0);
}

void Task256::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(x1, x2, x0, a1), 0.0);
  {
    // tensor label: I356
    std::unique_ptr<double[]> i0data = in(0)->get_block(a1, x0, x2, x1);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, a1.size(), x0.size(), x2.size(), x1.size());
  }
  out()->add_block(odata, x1, x2, x0, a1);
}

void Task257::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I356
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma52
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a1, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a1, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,2>(i1data, i1data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task258::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I356
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma49
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x0, x2, x1)]);
        sort_indices<0,1,2,3,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x3, a1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x3, a1)]);
        sort_indices<0,1,2,3,0,1,1,2>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), a1.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task259::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I356
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: Gamma50
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size());
    dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, x0.size()*x2.size()*x1.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task260::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a4, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c3, a4, c1, a2), 0.0);
  {
    // tensor label: I360
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
    sort_indices<2,1,0,3,1,1,1,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
  }
  out()->add_block(odata, c3, a4, c1, a2);
}

void Task261::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I360
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a4, c3, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a4, c3, a2), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<0,3,2,1,1,1,4,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->add_block(odata, c1, a4, c3, a2);
}

void Task262::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  {
    // tensor label: I362
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, c2, a1, x0);
    sort_indices<1,0,3,2,1,1,1,1>(i0data, odata, a3.size(), c2.size(), a1.size(), x0.size());
  }
  out()->add_block(odata, c2, a3, x0, a1);
}

void Task263::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I362
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, a1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, a1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, a1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, a1, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma29
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I363
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a3, c2, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), a3.size(), c2.size(), a1.size());
    dgemm_("T", "N", x0.size(), a3.size()*c2.size()*a1.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x0.size(), a3.size(), c2.size(), a1.size());
  out()->add_block(odata, a3, c2, a1, x0);
}

void Task264::Task_local::compute() {
  const Index x1 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I363
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, a3, c2, a1)]);
  std::fill_n(odata.get(), out()->get_size(x1, a3, c2, a1), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x1.size(), a3.size(), c2.size(), a1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x1, a1, c2, a3);
    sort_indices<0,3,2,1,1,1,2,1>(i1data, odata, x1.size(), a1.size(), c2.size(), a3.size());
  }
  out()->add_block(odata, x1, a3, c2, a1);
}

void Task265::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, a2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(x1, a2, x0, a1), 0.0);
  {
    // tensor label: I366
    std::unique_ptr<double[]> i0data = in(0)->get_block(a1, a2, x0, x1);
    sort_indices<3,1,2,0,1,1,1,1>(i0data, odata, a1.size(), a2.size(), x0.size(), x1.size());
  }
  out()->add_block(odata, x1, a2, x0, a1);
}

void Task266::Task_local::compute() {
  const Index a1 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I366
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, a2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, a2, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, a2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, a2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma50
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, x2, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), x2.size(), a2.size());
      dgemm_("T", "N", x0.size()*x1.size(), a1.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), a2.size());
  out()->add_block(odata, a1, a2, x0, x1);
}

#endif
