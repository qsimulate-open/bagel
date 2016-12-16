//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_tasks13.cc
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

#include <src/smith/mrci/MRCI_tasks13.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task600::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x7 = b(3);
  const Index a1 = b(4);
  const Index x6 = b(5);
  // tensor label: I1107
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x4, x3, x7, a1, x6)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x3, x7, a1, x6), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x5, x4, x3, x7, a1, x6)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x5, x4, x3, x7, a1, x6), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a1, x6, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a1, x6, a2)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x7.size(), a1.size(), x6.size(), a2.size());
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, a2, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, a2, x3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), a2.size(), x3.size());
    dgemm_("T", "N", x7.size()*a1.size()*x6.size(), x5.size()*x4.size()*x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x7.size()*a1.size()*x6.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x7.size(), a1.size(), x6.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, x5, x4, x3, x7, a1, x6);
}

void Task601::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: Gamma556
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x5, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x0, x6, x5, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x6.size(), x5.size(), x2.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x7, a1, x6, x5);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x7, a1, x6, x5)]);
        sort_indices<0,2,3,1,0,1,-1,1>(i1data, i1data_sorted, x7.size(), a1.size(), x6.size(), x5.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task602::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x9 : *range_[1]) {
    for (auto& x8 : *range_[1]) {
      for (auto& x7 : *range_[1]) {
        // tensor label: Gamma557
        std::unique_ptr<double[]> i0data = in(0)->get_block(x9, x0, x8, x7, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x9, x0, x8, x7, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x9.size(), x0.size(), x8.size(), x7.size(), x2.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x9, a1, x8, x7);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x9, a1, x8, x7)]);
        sort_indices<0,2,3,1,0,1,-1,2>(i1data, i1data_sorted, x9.size(), a1.size(), x8.size(), x7.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x9.size()*x8.size()*x7.size(),
               1.0, i0data_sorted, x9.size()*x8.size()*x7.size(), i1data_sorted, x9.size()*x8.size()*x7.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task603::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a4, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c3, a4, c1, a2), 0.0);
  {
    // tensor label: I162
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, c1, a4, c3);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, a2.size(), c1.size(), a4.size(), c3.size());
  }
  {
    // tensor label: I162
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3, a2, c1);
    sort_indices<1,0,3,2,1,1,1,1>(i0data, odata, a4.size(), c3.size(), a2.size(), c1.size());
  }
  out()->add_block(odata, c3, a4, c1, a2);
}

void Task604::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x1.size());
    // tensor label: I163
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a2.size(), x1.size());
    dgemm_("T", "N", c1.size()*a4.size()*c3.size(), a2.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size()*a4.size()*c3.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task605::Task_local::compute() {
  const Index a2 = b(0);
  const Index x1 = b(1);
  // tensor label: I163
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a2)]);
    sort_indices<0,1,0,1,-1,1>(i1data, i1data_sorted, x0.size(), a2.size());
    dgemm_("T", "N", x1.size(), a2.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), a2.size());
  out()->add_block(odata, a2, x1);
}

void Task606::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x1.size());
    // tensor label: I166
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a4.size(), x1.size());
    dgemm_("T", "N", c1.size()*a2.size()*c3.size(), a4.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task607::Task_local::compute() {
  const Index a4 = b(0);
  const Index x1 = b(1);
  // tensor label: I166
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, x1)]);
  std::fill_n(odata.get(), out()->get_size(a4, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4)]);
    sort_indices<0,1,0,1,2,1>(i1data, i1data_sorted, x0.size(), a4.size());
    dgemm_("T", "N", x1.size(), a4.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), a4.size());
  out()->add_block(odata, a4, x1);
}

void Task608::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  // tensor label: h1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size());
  // tensor label: I169
  std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c1)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), c1.size());
  dgemm_("T", "N", c3.size()*a2.size(), a4.size()*c1.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c3.size()*a2.size());
  sort_indices<1,3,2,0,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), a4.size(), c1.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task609::Task_local::compute() {
  const Index a4 = b(0);
  const Index c1 = b(1);
  // tensor label: I169
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, c1)]);
  std::fill_n(odata.get(), out()->get_size(a4, c1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma32
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: I170
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a4, c1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a4, c1, x0)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), a4.size(), c1.size(), x0.size());
      dgemm_("T", "N", 1, a4.size()*c1.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size());
  out()->add_block(odata, a4, c1);
}

void Task610::Task_local::compute() {
  const Index x1 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: I170
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, a4, c1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, a4, c1, x0), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c1, x0);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x1.size(), a4.size(), c1.size(), x0.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a4, x1, x0);
    sort_indices<2,1,0,3,1,1,-2,1>(i1data, odata, c1.size(), a4.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x1, a4, c1, x0);
}

void Task611::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  // tensor label: h1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size());
  // tensor label: I172
  std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c1)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a2.size(), c1.size());
  dgemm_("T", "N", c3.size()*a4.size(), a2.size()*c1.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c3.size()*a4.size());
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, c3.size(), a4.size(), a2.size(), c1.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task612::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  // tensor label: I172
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma32
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: I173
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a2, c1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a2, c1, x0)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), a2.size(), c1.size(), x0.size());
      dgemm_("T", "N", 1, a2.size()*c1.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size());
  out()->add_block(odata, a2, c1);
}

void Task613::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: I173
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, a2, c1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, a2, c1, x0), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c1, x0);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x1.size(), a2.size(), c1.size(), x0.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, x1, x0);
    sort_indices<2,1,0,3,1,1,4,1>(i1data, odata, c1.size(), a2.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x1, a2, c1, x0);
}

void Task614::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& c5 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c5, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c5, a2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c5.size(), a2.size());
    // tensor label: I181
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, c5);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, c5)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), c5.size());
    dgemm_("T", "N", c1.size()*a4.size()*a2.size(), c3.size(), c5.size(),
           1.0, i0data_sorted, c5.size(), i1data_sorted, c5.size(),
           1.0, odata_sorted, c1.size()*a4.size()*a2.size());
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), a2.size(), c3.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task615::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I181
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, c5)]);
  std::fill_n(odata.get(), out()->get_size(c3, c5), 0.0);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c5);
    sort_indices<0,1,1,1,4,1>(i0data, odata, c3.size(), c5.size());
  }
  out()->add_block(odata, c3, c5);
}

void Task616::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I181
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, c5)]);
  std::fill_n(odata.get(), out()->get_size(c3, c5), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, c5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c5), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma32
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: I1243
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, c5, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, c5, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), c5.size(), x1.size(), x0.size());
      dgemm_("T", "N", 1, c3.size()*c5.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c3.size(), c5.size());
  out()->add_block(odata, c3, c5);
}

void Task617::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1243
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, c5, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, c5, x1, x0), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c5, x1, x0);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, c3.size(), c5.size(), x1.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x1, c5, c3, x0);
    sort_indices<2,1,0,3,1,1,-1,1>(i1data, odata, x1.size(), c5.size(), c3.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i2data = in(0)->get_block(x1, x0, c3, c5);
    sort_indices<2,3,0,1,1,1,2,1>(i2data, odata, x1.size(), x0.size(), c3.size(), c5.size());
  }
  out()->add_block(odata, c3, c5, x1, x0);
}

void Task618::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I181
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, c5)]);
  std::fill_n(odata.get(), out()->get_size(c3, c5), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, c5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c5), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma12
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x1, x0, c5);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x1, x0, c5)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), x1.size(), x0.size(), c5.size());
      dgemm_("T", "N", 1, c3.size()*c5.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c3.size(), c5.size());
  out()->add_block(odata, c3, c5);
}

void Task619::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& c5 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c5, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c5, a4)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c5.size(), a4.size());
    // tensor label: I183
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, c5);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, c5)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), c5.size());
    dgemm_("T", "N", c1.size()*a2.size()*a4.size(), c3.size(), c5.size(),
           1.0, i0data_sorted, c5.size(), i1data_sorted, c5.size(),
           1.0, odata_sorted, c1.size()*a2.size()*a4.size());
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), a4.size(), c3.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task620::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I183
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, c5)]);
  std::fill_n(odata.get(), out()->get_size(c3, c5), 0.0);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c5);
    sort_indices<0,1,1,1,-8,1>(i0data, odata, c3.size(), c5.size());
  }
  out()->add_block(odata, c3, c5);
}

void Task621::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I183
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, c5)]);
  std::fill_n(odata.get(), out()->get_size(c3, c5), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, c5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c5), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma32
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: I1246
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, c5, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, c5, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), c5.size(), x1.size(), x0.size());
      dgemm_("T", "N", 1, c3.size()*c5.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c3.size(), c5.size());
  out()->add_block(odata, c3, c5);
}

void Task622::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1246
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, c5, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, c5, x1, x0), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c5, x1, x0);
    sort_indices<0,1,2,3,1,1,-4,1>(i0data, odata, c3.size(), c5.size(), x1.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x1, c5, c3, x0);
    sort_indices<2,1,0,3,1,1,2,1>(i1data, odata, x1.size(), c5.size(), c3.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i2data = in(0)->get_block(x1, x0, c3, c5);
    sort_indices<2,3,0,1,1,1,-4,1>(i2data, odata, x1.size(), x0.size(), c3.size(), c5.size());
  }
  out()->add_block(odata, c3, c5, x1, x0);
}

void Task623::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I183
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, c5)]);
  std::fill_n(odata.get(), out()->get_size(c3, c5), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, c5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c5), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma12
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x1, x0, c5);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x1, x0, c5)]);
      sort_indices<2,1,0,3,0,1,-2,1>(i1data, i1data_sorted, c3.size(), x1.size(), x0.size(), c5.size());
      dgemm_("T", "N", 1, c3.size()*c5.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c3.size(), c5.size());
  out()->add_block(odata, c3, c5);
}

void Task624::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& a5 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a5, c3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a5, c3, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a5.size(), c3.size(), a2.size());
    // tensor label: I185
    std::unique_ptr<double[]> i1data = in(1)->get_block(a5, a4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, a4)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a5.size(), a4.size());
    dgemm_("T", "N", c1.size()*c3.size()*a2.size(), a4.size(), a5.size(),
           1.0, i0data_sorted, a5.size(), i1data_sorted, a5.size(),
           1.0, odata_sorted, c1.size()*c3.size()*a2.size());
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), a2.size(), a4.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task625::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  // tensor label: I185
  std::unique_ptr<double[]> odata(new double[out()->get_size(a5, a4)]);
  std::fill_n(odata.get(), out()->get_size(a5, a4), 0.0);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, a4);
    sort_indices<0,1,1,1,-4,1>(i0data, odata, a5.size(), a4.size());
  }
  out()->add_block(odata, a5, a4);
}

void Task626::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  // tensor label: I185
  std::unique_ptr<double[]> odata(new double[out()->get_size(a5, a4)]);
  std::fill_n(odata.get(), out()->get_size(a5, a4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a5, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a5, a4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma32
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: I1249
      std::unique_ptr<double[]> i1data = in(1)->get_block(a5, a4, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, a4, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a5.size(), a4.size(), x1.size(), x0.size());
      dgemm_("T", "N", 1, a5.size()*a4.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a5.size(), a4.size());
  out()->add_block(odata, a5, a4);
}

void Task627::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1249
  std::unique_ptr<double[]> odata(new double[out()->get_size(a5, a4, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a5, a4, x1, x0), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, a4, x1, x0);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, a5.size(), a4.size(), x1.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x1, a4, a5, x0);
    sort_indices<2,1,0,3,1,1,1,1>(i1data, odata, x1.size(), a4.size(), a5.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i2data = in(0)->get_block(x1, x0, a5, a4);
    sort_indices<2,3,0,1,1,1,-2,1>(i2data, odata, x1.size(), x0.size(), a5.size(), a4.size());
  }
  out()->add_block(odata, a5, a4, x1, x0);
}

void Task628::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  // tensor label: I185
  std::unique_ptr<double[]> odata(new double[out()->get_size(a5, a4)]);
  std::fill_n(odata.get(), out()->get_size(a5, a4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a5, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a5, a4), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma12
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(a5, x1, x0, a4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, x1, x0, a4)]);
      sort_indices<2,1,0,3,0,1,-1,1>(i1data, i1data_sorted, a5.size(), x1.size(), x0.size(), a4.size());
      dgemm_("T", "N", 1, a5.size()*a4.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a5.size(), a4.size());
  out()->add_block(odata, a5, a4);
}

void Task629::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& a5 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a5);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a5)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a5.size());
    // tensor label: I187
    std::unique_ptr<double[]> i1data = in(1)->get_block(a5, a4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, a4)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a5.size(), a4.size());
    dgemm_("T", "N", c1.size()*a2.size()*c3.size(), a4.size(), a5.size(),
           1.0, i0data_sorted, a5.size(), i1data_sorted, a5.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task630::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  // tensor label: I187
  std::unique_ptr<double[]> odata(new double[out()->get_size(a5, a4)]);
  std::fill_n(odata.get(), out()->get_size(a5, a4), 0.0);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, a4);
    sort_indices<0,1,1,1,8,1>(i0data, odata, a5.size(), a4.size());
  }
  out()->add_block(odata, a5, a4);
}

void Task631::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  // tensor label: I187
  std::unique_ptr<double[]> odata(new double[out()->get_size(a5, a4)]);
  std::fill_n(odata.get(), out()->get_size(a5, a4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a5, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a5, a4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma32
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: I1252
      std::unique_ptr<double[]> i1data = in(1)->get_block(a5, a4, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, a4, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a5.size(), a4.size(), x1.size(), x0.size());
      dgemm_("T", "N", 1, a5.size()*a4.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a5.size(), a4.size());
  out()->add_block(odata, a5, a4);
}

void Task632::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1252
  std::unique_ptr<double[]> odata(new double[out()->get_size(a5, a4, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a5, a4, x1, x0), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, a4, x1, x0);
    sort_indices<0,1,2,3,1,1,4,1>(i0data, odata, a5.size(), a4.size(), x1.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x1, a4, a5, x0);
    sort_indices<2,1,0,3,1,1,-2,1>(i1data, odata, x1.size(), a4.size(), a5.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i2data = in(0)->get_block(x1, x0, a5, a4);
    sort_indices<2,3,0,1,1,1,4,1>(i2data, odata, x1.size(), x0.size(), a5.size(), a4.size());
  }
  out()->add_block(odata, a5, a4, x1, x0);
}

void Task633::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  // tensor label: I187
  std::unique_ptr<double[]> odata(new double[out()->get_size(a5, a4)]);
  std::fill_n(odata.get(), out()->get_size(a5, a4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a5, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a5, a4), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma12
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(a5, x1, x0, a4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, x1, x0, a4)]);
      sort_indices<2,1,0,3,0,1,2,1>(i1data, i1data_sorted, a5.size(), x1.size(), x0.size(), a4.size());
      dgemm_("T", "N", 1, a5.size()*a4.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a5.size(), a4.size());
  out()->add_block(odata, a5, a4);
}

void Task634::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c1, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c1.size(), a2.size());
    // tensor label: I189
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), x1.size());
    dgemm_("T", "N", a4.size()*c1.size()*a2.size(), c3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a4.size()*c1.size()*a2.size());
  }
  sort_indices<2,1,0,3,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), a2.size(), c3.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task635::Task_local::compute() {
  const Index c3 = b(0);
  const Index x1 = b(1);
  // tensor label: I189
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c3, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma32
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x0)]);
    sort_indices<1,0,0,1,-2,1>(i1data, i1data_sorted, c3.size(), x0.size());
    dgemm_("T", "N", x1.size(), c3.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), c3.size());
  out()->add_block(odata, c3, x1);
}

void Task636::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c1, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2, c1, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), a4.size());
    // tensor label: I192
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), x1.size());
    dgemm_("T", "N", a2.size()*c1.size()*a4.size(), c3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a2.size()*c1.size()*a4.size());
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), a4.size(), c3.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task637::Task_local::compute() {
  const Index c3 = b(0);
  const Index x1 = b(1);
  // tensor label: I192
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c3, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma32
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x0)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), x0.size());
    dgemm_("T", "N", x1.size(), c3.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), c3.size());
  out()->add_block(odata, c3, x1);
}

void Task638::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a4, c3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a4, c3, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a4.size(), c3.size(), a2.size());
    // tensor label: I1116
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x0)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), x0.size());
    dgemm_("T", "N", a4.size()*c3.size()*a2.size(), c1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a4.size()*c3.size()*a2.size());
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a4.size(), c3.size(), a2.size(), c1.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task639::Task_local::compute() {
  const Index c1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1116
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma10
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, c1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, c1, x1)]);
        sort_indices<0,1,3,2,0,1,-1,1>(i1data, i1data_sorted, x3.size(), x2.size(), c1.size(), x1.size());
        dgemm_("T", "N", x0.size(), c1.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c1.size());
  out()->add_block(odata, c1, x0);
}

void Task640::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a2, c3, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a2, c3, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a2.size(), c3.size(), a4.size());
    // tensor label: I1119
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x0)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), x0.size());
    dgemm_("T", "N", a2.size()*c3.size()*a4.size(), c1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a2.size()*c3.size()*a4.size());
  }
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, a2.size(), c3.size(), a4.size(), c1.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task641::Task_local::compute() {
  const Index c1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1119
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma10
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, c1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, c1, x1)]);
        sort_indices<0,1,3,2,0,1,2,1>(i1data, i1data_sorted, x3.size(), x2.size(), c1.size(), x1.size());
        dgemm_("T", "N", x0.size(), c1.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c1.size());
  out()->add_block(odata, c1, x0);
}

void Task642::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x3.size());
    // tensor label: I1122
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a2.size(), x3.size());
    dgemm_("T", "N", c1.size()*a4.size()*c3.size(), a2.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c1.size()*a4.size()*c3.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task643::Task_local::compute() {
  const Index a2 = b(0);
  const Index x3 = b(1);
  // tensor label: I1122
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, x3)]);
  std::fill_n(odata.get(), out()->get_size(a2, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: Gamma5
        std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x3, x1, x0)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a2, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a2, x1, x0)]);
        sort_indices<0,2,3,1,0,1,-1,2>(i1data, i1data_sorted, x2.size(), a2.size(), x1.size(), x0.size());
        dgemm_("T", "N", x3.size(), a2.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), a2.size());
  out()->add_block(odata, a2, x3);
}

void Task644::Task_local::compute() {
  const Index a2 = b(0);
  const Index x3 = b(1);
  // tensor label: I1122
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, x3)]);
  std::fill_n(odata.get(), out()->get_size(a2, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x3), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma197
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x3, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, x0, a2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, x0, a2)]);
        sort_indices<2,0,1,3,0,1,-1,2>(i1data, i1data_sorted, x2.size(), x1.size(), x0.size(), a2.size());
        dgemm_("T", "N", x3.size(), a2.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), a2.size());
  out()->add_block(odata, a2, x3);
}

void Task645::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x3.size());
    // tensor label: I1125
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a4.size(), x3.size());
    dgemm_("T", "N", c1.size()*a2.size()*c3.size(), a4.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task646::Task_local::compute() {
  const Index a4 = b(0);
  const Index x3 = b(1);
  // tensor label: I1125
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, x3)]);
  std::fill_n(odata.get(), out()->get_size(a4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: Gamma5
        std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x3, x1, x0)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a4, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a4, x1, x0)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x2.size(), a4.size(), x1.size(), x0.size());
        dgemm_("T", "N", x3.size(), a4.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), a4.size());
  out()->add_block(odata, a4, x3);
}

void Task647::Task_local::compute() {
  const Index a4 = b(0);
  const Index x3 = b(1);
  // tensor label: I1125
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, x3)]);
  std::fill_n(odata.get(), out()->get_size(a4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x3), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma197
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x3, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, x0, a4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, x0, a4)]);
        sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), x0.size(), a4.size());
        dgemm_("T", "N", x3.size(), a4.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), a4.size());
  out()->add_block(odata, a4, x3);
}

void Task648::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& c5 : *range_[0]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c5, a4, c1, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c5, a4, c1, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, c5.size(), a4.size(), c1.size(), x1.size());
      // tensor label: I1134
      std::unique_ptr<double[]> i1data = in(1)->get_block(c5, c3, a2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c5, c3, a2, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, c5.size(), c3.size(), a2.size(), x1.size());
      dgemm_("T", "N", a4.size()*c1.size(), c3.size()*a2.size(), c5.size()*x1.size(),
             1.0, i0data_sorted, c5.size()*x1.size(), i1data_sorted, c5.size()*x1.size(),
             1.0, odata_sorted, a4.size()*c1.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), c3.size(), a2.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task649::Task_local::compute() {
  const Index c5 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1134
  std::unique_ptr<double[]> odata(new double[out()->get_size(c5, c3, a2, x1)]);
  std::fill_n(odata.get(), out()->get_size(c5, c3, a2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c5, c3, a2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c5, c3, a2, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I1135
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c5, c3, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c5, c3, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c5.size(), c3.size(), a2.size());
    dgemm_("T", "N", x1.size(), c5.size()*c3.size()*a2.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), c5.size(), c3.size(), a2.size());
  out()->add_block(odata, c5, c3, a2, x1);
}

#endif
