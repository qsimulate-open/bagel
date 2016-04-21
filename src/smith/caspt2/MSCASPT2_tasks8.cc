//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_tasks8.cc
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

#include <src/smith/caspt2/MSCASPT2_tasks8.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MSCASPT2;

void Task350::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index c1 = b(3);
  // tensor label: I348
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x4, x3, c1)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x3, c1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x5, x4, x3, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x5, x4, x3, c1), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, c2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), c2.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, c2, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, c2, x3)]);
    sort_indices<2,0,1,3,0,1,-1,1>(i1data, i1data_sorted, x5.size(), x4.size(), c2.size(), x3.size());
    dgemm_("T", "N", c1.size(), x5.size()*x4.size()*x3.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, c1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c1.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, x5, x4, x3, c1);
}

void Task351::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index c1 = b(3);
  // tensor label: I348
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x4, x3, c1)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x3, c1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x5, x4, x3, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x5, x4, x3, c1), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a2, x3)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a2.size(), x3.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, x5, x4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, x5, x4)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), x5.size(), x4.size());
    dgemm_("T", "N", x3.size(), c1.size()*x5.size()*x4.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x3.size());
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), c1.size(), x5.size(), x4.size());
  out()->add_block(odata, x5, x4, x3, c1);
}

void Task352::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I347
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x4, x3, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x3, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x5, x4, x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x5, x4, x3, x2, x1, x0), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c1, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, c1, x3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), c1.size(), x3.size());
    // tensor label: I488
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, c1, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, c1, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), c1.size(), x2.size());
    dgemm_("T", "N", x5.size()*x4.size()*x3.size(), x1.size()*x0.size()*x2.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x5.size()*x4.size()*x3.size());
  }
  sort_indices<0,1,2,5,3,4,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), x1.size(), x0.size(), x2.size());
  out()->add_block(odata, x5, x4, x3, x2, x1, x0);
}

void Task353::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index x2 = b(3);
  // tensor label: I488
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0, c1, x2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, c1, x2), 0.0);
  {
    // tensor label: l2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c1, x2);
    dscal_(x2.size()*c1.size()*x1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<1,0,2,3,1,1,-1,1>(i0data, odata, x0.size(), x1.size(), c1.size(), x2.size());
  }
  out()->add_block(odata, x1, x0, c1, x2);
}

void Task354::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index x2 = b(3);
  // tensor label: I488
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0, c1, x2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, c1, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, x2), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x2.size(), a2.size());
    // tensor label: l2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, x0, x1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), x0.size(), x1.size());
    dgemm_("T", "N", x2.size(), x1.size()*x0.size()*c1.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x2.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x2.size(), c1.size(), x0.size(), x1.size());
  out()->add_block(odata, x1, x0, c1, x2);
}

void Task355::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I347
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x4, x3, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x3, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x5, x4, x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x5, x4, x3, x2, x1, x0), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: l2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c1, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, c1, x3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), c1.size(), x3.size());
    // tensor label: I710
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, x2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, x2, c1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), x2.size(), c1.size());
    dgemm_("T", "N", x3.size()*x4.size()*x5.size(), x0.size()*x1.size()*x2.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x3.size()*x4.size()*x5.size());
  }
  sort_indices<0,1,2,5,4,3,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), x0.size(), x1.size(), x2.size());
  out()->add_block(odata, x5, x4, x3, x2, x1, x0);
}

void Task356::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index c1 = b(3);
  // tensor label: I710
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1, x2, c1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1, x2, c1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, x2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, x2, c1), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, c2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), c2.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, c2, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, c2, x2)]);
    sort_indices<2,0,1,3,0,1,-1,1>(i1data, i1data_sorted, x0.size(), x1.size(), c2.size(), x2.size());
    dgemm_("T", "N", c1.size(), x0.size()*x1.size()*x2.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, c1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c1.size(), x0.size(), x1.size(), x2.size());
  out()->add_block(odata, x0, x1, x2, c1);
}

void Task357::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index c1 = b(3);
  // tensor label: I710
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1, x2, c1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1, x2, c1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, x2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, x2, c1), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a2, x2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a2.size(), x2.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, x0, x1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), x0.size(), x1.size());
    dgemm_("T", "N", x2.size(), c1.size()*x0.size()*x1.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x2.size());
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x2.size(), c1.size(), x0.size(), x1.size());
  out()->add_block(odata, x0, x1, x2, c1);
}

void Task358::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I347
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x4, x3, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x3, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x5, x4, x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x5, x4, x3, x2, x1, x0), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1, c1, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size(), c1.size(), x2.size());
    // tensor label: I850
    std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x5, c1, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x5, c1, x3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x4.size(), x5.size(), c1.size(), x3.size());
    dgemm_("T", "N", x0.size()*x1.size()*x2.size(), x4.size()*x5.size()*x3.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x0.size()*x1.size()*x2.size());
  }
  sort_indices<4,3,5,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x4.size(), x5.size(), x3.size());
  out()->add_block(odata, x5, x4, x3, x2, x1, x0);
}

void Task359::Task_local::compute() {
  const Index x4 = b(0);
  const Index x5 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I850
  std::unique_ptr<double[]> odata(new double[out()->get_size(x4, x5, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(x4, x5, c1, x3), 0.0);
  {
    // tensor label: l2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c1, x3);
    dscal_(x3.size()*c1.size()*x4.size()*x5.size(), e0_, i0data.get(), 1);
    sort_indices<1,0,2,3,1,1,-1,1>(i0data, odata, x5.size(), x4.size(), c1.size(), x3.size());
  }
  out()->add_block(odata, x4, x5, c1, x3);
}

void Task360::Task_local::compute() {
  const Index x4 = b(0);
  const Index x5 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I850
  std::unique_ptr<double[]> odata(new double[out()->get_size(x4, x5, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(x4, x5, c1, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, x5, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x5, c1, x3), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size());
    // tensor label: l2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, x5, x4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, x5, x4)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), x5.size(), x4.size());
    dgemm_("T", "N", x3.size(), x4.size()*x5.size()*c1.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x3.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x3.size(), c1.size(), x5.size(), x4.size());
  out()->add_block(odata, x4, x5, c1, x3);
}

void Task361::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I324
  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        for (auto& x0 : *range_[1]) {
          // tensor label: Gamma117
          std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x2, x3, x1, x0);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x2, x3, x1, x0)]);
          sort_indices<1,2,3,4,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x2.size(), x3.size(), x1.size(), x0.size());
          // tensor label: I351
          std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
          sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
          dgemm_("T", "N", ci0.size(), 1, x3.size()*x2.size()*x1.size()*x0.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x1.size()*x0.size(), i1data_sorted, x3.size()*x2.size()*x1.size()*x0.size(),
                 1.0, odata_sorted, ci0.size());
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->add_block(odata, ci0);
}

void Task362::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I351
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: l2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1, c1, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size(), c1.size(), x2.size());
    // tensor label: I352
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x3)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), x3.size());
    dgemm_("T", "N", x2.size()*x1.size()*x0.size(), x3.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x2.size()*x1.size()*x0.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size());
  out()->add_block(odata, x3, x2, x1, x0);
}

void Task363::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  // tensor label: I352
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(c1, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(a3, c2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, c2)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), c2.size());
      // tensor label: I353
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a3, c1, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a3, c1, x3)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), a3.size(), c1.size(), x3.size());
      dgemm_("T", "N", 1, c1.size()*x3.size(), c2.size()*a3.size(),
             1.0, i0data_sorted, c2.size()*a3.size(), i1data_sorted, c2.size()*a3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size());
  out()->add_block(odata, c1, x3);
}

void Task364::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I353
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, a3, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(c2, a3, c1, x3), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c1, x3);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, c2.size(), a3.size(), c1.size(), x3.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a3, c2, x3);
    sort_indices<2,1,0,3,1,1,-1,1>(i1data, odata, c1.size(), a3.size(), c2.size(), x3.size());
  }
  out()->add_block(odata, c2, a3, c1, x3);
}

void Task365::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I351
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: l2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, x1)]);
      sort_indices<2,1,0,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), x1.size());
      // tensor label: I442
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, x3, x2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x3.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), x3.size()*x2.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x3.size(), x2.size());
  out()->add_block(odata, x3, x2, x1, x0);
}

void Task366::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I442
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, a1, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c2, a1, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x3, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, c3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x2.size(), c3.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, c3, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, c3, x3)]);
    sort_indices<2,0,1,3,0,1,-1,1>(i1data, i1data_sorted, c2.size(), a1.size(), c3.size(), x3.size());
    dgemm_("T", "N", x2.size(), c2.size()*a1.size()*x3.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, x2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x2.size(), c2.size(), a1.size(), x3.size());
  out()->add_block(odata, c2, a1, x3, x2);
}

void Task367::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I351
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      // tensor label: l2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x0, x1)]);
      sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x0.size(), x1.size());
      // tensor label: I492
      std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c1, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c1, x3, x2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), c1.size(), x3.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), x3.size()*x2.size(), a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x3.size(), x2.size());
  out()->add_block(odata, x3, x2, x1, x0);
}

void Task368::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I492
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, x3, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, c3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x2.size(), c3.size());
    // tensor label: I493
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, c1, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, c1, x3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), c1.size(), x3.size());
    dgemm_("T", "N", x2.size(), a2.size()*c1.size()*x3.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, x2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x2.size(), a2.size(), c1.size(), x3.size());
  out()->add_block(odata, a2, c1, x3, x2);
}

void Task369::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I493
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a2, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, c1, x3), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x3);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c3.size(), a2.size(), c1.size(), x3.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, x3);
    sort_indices<2,1,0,3,1,1,2,1>(i1data, odata, c1.size(), a2.size(), c3.size(), x3.size());
  }
  out()->add_block(odata, c3, a2, c1, x3);
}

void Task370::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I351
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c3, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1, c3, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size(), c3.size(), x2.size());
    // tensor label: I734
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x3.size(), c3.size());
    dgemm_("T", "N", x0.size()*x1.size()*x2.size(), x3.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, x0.size()*x1.size()*x2.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size());
  out()->add_block(odata, x3, x2, x1, x0);
}

void Task371::Task_local::compute() {
  const Index x3 = b(0);
  const Index c3 = b(1);
  // tensor label: I734
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, c3)]);
  std::fill_n(odata.get(), out()->get_size(x3, c3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, c3), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size());
      // tensor label: l2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x3)]);
      sort_indices<0,1,2,3,0,1,2,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x3.size());
      dgemm_("T", "N", 1, x3.size()*c3.size(), a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c3.size(), x3.size());
  out()->add_block(odata, x3, c3);
}

void Task372::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I351
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1, c1, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size(), c1.size(), x2.size());
    // tensor label: I738
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x3.size(), c1.size());
    dgemm_("T", "N", x0.size()*x1.size()*x2.size(), x3.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x0.size()*x1.size()*x2.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size());
  out()->add_block(odata, x3, x2, x1, x0);
}

void Task373::Task_local::compute() {
  const Index x3 = b(0);
  const Index c1 = b(1);
  // tensor label: I738
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, c1)]);
  std::fill_n(odata.get(), out()->get_size(x3, c1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, c1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size());
      // tensor label: l2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x3)]);
      sort_indices<2,1,0,3,0,1,-1,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x3.size());
      dgemm_("T", "N", 1, x3.size()*c1.size(), c3.size()*a2.size(),
             1.0, i0data_sorted, c3.size()*a2.size(), i1data_sorted, c3.size()*a2.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size());
  out()->add_block(odata, x3, c1);
}

void Task374::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I351
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a2, c1, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a2, c1, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a2.size(), c1.size(), x1.size());
      // tensor label: I776
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a2, c1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a2, c1, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a2.size(), c1.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), x3.size()*x2.size(), a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x3.size(), x2.size());
  out()->add_block(odata, x3, x2, x1, x0);
}

void Task375::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x2 = b(3);
  // tensor label: I776
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a2, c1, x2)]);
  std::fill_n(odata.get(), out()->get_size(x3, a2, c1, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a2, c1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a2, c1, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), x2.size());
    // tensor label: l2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x3)]);
    sort_indices<2,0,1,3,0,1,-1,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x3.size());
    dgemm_("T", "N", x2.size(), x3.size()*a2.size()*c1.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, x2.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x2.size(), c1.size(), a2.size(), x3.size());
  out()->add_block(odata, x3, a2, c1, x2);
}

void Task376::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I351
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), x0.size(), x1.size());
      // tensor label: I780
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c3, a2, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c3, a2, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), c3.size(), a2.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), x3.size()*x2.size(), c3.size()*a2.size(),
             1.0, i0data_sorted, c3.size()*a2.size(), i1data_sorted, c3.size()*a2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x3.size(), x2.size());
  out()->add_block(odata, x3, x2, x1, x0);
}

void Task377::Task_local::compute() {
  const Index x3 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I780
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, c3, a2, x2)]);
  std::fill_n(odata.get(), out()->get_size(x3, c3, a2, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, c3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, c3, a2, x2), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), x2.size());
    // tensor label: l2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x3)]);
    sort_indices<0,1,2,3,0,1,-1,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x3.size());
    dgemm_("T", "N", x2.size(), x3.size()*c3.size()*a2.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x2.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x2.size(), a2.size(), c3.size(), x3.size());
  out()->add_block(odata, x3, c3, a2, x2);
}

void Task378::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I351
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x0.size(), x1.size());
      // tensor label: I784
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a2, c1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a2, c1, x2)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a2.size(), c1.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), x3.size()*x2.size(), a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x3.size(), x2.size());
  out()->add_block(odata, x3, x2, x1, x0);
}

void Task379::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x2 = b(3);
  // tensor label: I784
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a2, c1, x2)]);
  std::fill_n(odata.get(), out()->get_size(x3, a2, c1, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a2, c1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a2, c1, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), x2.size());
    // tensor label: l2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x3)]);
    sort_indices<2,0,1,3,0,1,2,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x3.size());
    dgemm_("T", "N", x2.size(), x3.size()*a2.size()*c1.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, x2.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x2.size(), c1.size(), a2.size(), x3.size());
  out()->add_block(odata, x3, a2, c1, x2);
}

void Task380::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I324
  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        for (auto& x4 : *range_[1]) {
          for (auto& x1 : *range_[1]) {
            for (auto& x0 : *range_[1]) {
              // tensor label: Gamma119
              std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x3, x2, x4, x1, x0);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x5, x3, x2, x4, x1, x0)]);
              sort_indices<1,2,3,4,5,6,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x5.size(), x3.size(), x2.size(), x4.size(), x1.size(), x0.size());
              // tensor label: I359
              std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x3, x2, x1, x0);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x3, x2, x1, x0)]);
              sort_indices<0,2,3,1,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());
              dgemm_("T", "N", ci0.size(), 1, x5.size()*x4.size()*x3.size()*x2.size()*x1.size()*x0.size(),
                     1.0, i0data_sorted, x5.size()*x4.size()*x3.size()*x2.size()*x1.size()*x0.size(), i1data_sorted, x5.size()*x4.size()*x3.size()*x2.size()*x1.size()*x0.size(),
                     1.0, odata_sorted, ci0.size());
            }
          }
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->add_block(odata, ci0);
}

void Task381::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I359
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x4, x3, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x3, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x5, x4, x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x5, x4, x3, x2, x1, x0), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: l2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1, c1, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size(), c1.size(), x2.size());
    // tensor label: I360
    std::unique_ptr<double[]> i1data = in(1)->get_block(x5, c1, x4, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, c1, x4, x3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), c1.size(), x4.size(), x3.size());
    dgemm_("T", "N", x2.size()*x1.size()*x0.size(), x5.size()*x4.size()*x3.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x2.size()*x1.size()*x0.size());
  }
  sort_indices<3,4,5,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, x5, x4, x3, x2, x1, x0);
}

void Task382::Task_local::compute() {
  const Index x5 = b(0);
  const Index c1 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I360
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, c1, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(x5, c1, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x5, c1, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x5, c1, x4, x3), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a2, x3)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a2.size(), x3.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a2, c1, x4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a2, c1, x4)]);
    sort_indices<1,0,2,3,0,1,-1,1>(i1data, i1data_sorted, x5.size(), a2.size(), c1.size(), x4.size());
    dgemm_("T", "N", x3.size(), x5.size()*c1.size()*x4.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x3.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x3.size(), x5.size(), c1.size(), x4.size());
  out()->add_block(odata, x5, c1, x4, x3);
}

void Task383::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I359
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x4, x3, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x3, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x5, x4, x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x5, x4, x3, x2, x1, x0), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c2, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1, c2, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size(), c2.size(), x2.size());
    // tensor label: I796
    std::unique_ptr<double[]> i1data = in(1)->get_block(x4, c2, x5, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, c2, x5, x3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x4.size(), c2.size(), x5.size(), x3.size());
    dgemm_("T", "N", x0.size()*x1.size()*x2.size(), x4.size()*x5.size()*x3.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, x0.size()*x1.size()*x2.size());
  }
  sort_indices<4,3,5,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x4.size(), x5.size(), x3.size());
  out()->add_block(odata, x5, x4, x3, x2, x1, x0);
}

void Task384::Task_local::compute() {
  const Index x4 = b(0);
  const Index c2 = b(1);
  const Index x5 = b(2);
  const Index x3 = b(3);
  // tensor label: I796
  std::unique_ptr<double[]> odata(new double[out()->get_size(x4, c2, x5, x3)]);
  std::fill_n(odata.get(), out()->get_size(x4, c2, x5, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, c2, x5, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, c2, x5, x3), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size());
    // tensor label: l2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a1, c2, x4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a1, c2, x4)]);
    sort_indices<1,0,2,3,0,1,-1,1>(i1data, i1data_sorted, x5.size(), a1.size(), c2.size(), x4.size());
    dgemm_("T", "N", x3.size(), x4.size()*c2.size()*x5.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x3.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x3.size(), x5.size(), c2.size(), x4.size());
  out()->add_block(odata, x4, c2, x5, x3);
}

void Task385::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I324
  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        for (auto& x1 : *range_[1]) {
          // tensor label: Gamma122
          std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x0, x1);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x3, x2, x0, x1)]);
          sort_indices<1,2,3,4,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x3.size(), x2.size(), x0.size(), x1.size());
          // tensor label: I371
          std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x3, x2, x1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x3, x2, x1)]);
          sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x3.size(), x2.size(), x1.size());
          dgemm_("T", "N", ci0.size(), 1, x0.size()*x3.size()*x2.size()*x1.size(),
                 1.0, i0data_sorted, x0.size()*x3.size()*x2.size()*x1.size(), i1data_sorted, x0.size()*x3.size()*x2.size()*x1.size(),
                 1.0, odata_sorted, ci0.size());
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->add_block(odata, ci0);
}

void Task386::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I371
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x3, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x3, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3, x2, x1), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, c3, x1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), c3.size(), x1.size());
    // tensor label: I372
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size());
    dgemm_("T", "N", x3.size()*x2.size()*x1.size(), x0.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, x3.size()*x2.size()*x1.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size());
  out()->add_block(odata, x0, x3, x2, x1);
}

void Task387::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  // tensor label: I372
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, c3)]);
  std::fill_n(odata.get(), out()->get_size(x0, c3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size());
      // tensor label: l2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x0)]);
      sort_indices<0,1,2,3,0,1,2,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
      dgemm_("T", "N", 1, x0.size()*c3.size(), a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c3.size(), x0.size());
  out()->add_block(odata, x0, c3);
}

void Task388::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I371
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x3, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x3, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3, x2, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c1, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, c1, x1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), c1.size(), x1.size());
    // tensor label: I376
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), c1.size());
    dgemm_("T", "N", x3.size()*x2.size()*x1.size(), x0.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x3.size()*x2.size()*x1.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size());
  out()->add_block(odata, x0, x3, x2, x1);
}

void Task389::Task_local::compute() {
  const Index x0 = b(0);
  const Index c1 = b(1);
  // tensor label: I376
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, c1)]);
  std::fill_n(odata.get(), out()->get_size(x0, c1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size());
      // tensor label: l2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x0)]);
      sort_indices<2,1,0,3,0,1,-1,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
      dgemm_("T", "N", 1, x0.size()*c1.size(), c3.size()*a2.size(),
             1.0, i0data_sorted, c3.size()*a2.size(), i1data_sorted, c3.size()*a2.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c1.size(), x0.size());
  out()->add_block(odata, x0, c1);
}

void Task390::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I371
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x3, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x3, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3, x2, x1), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c1, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c1.size(), x2.size());
      // tensor label: I414
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a2, c1, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a2, c1, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a2.size(), c1.size(), x1.size());
      dgemm_("T", "N", x3.size()*x2.size(), x0.size()*x1.size(), a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x0.size(), x1.size());
  out()->add_block(odata, x0, x3, x2, x1);
}

void Task391::Task_local::compute() {
  const Index x0 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x1 = b(3);
  // tensor label: I414
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, a2, c1, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, a2, c1, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a2, c1, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a2, c1, x1), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), x1.size());
    // tensor label: l2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x0)]);
    sort_indices<2,0,1,3,0,1,-1,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
    dgemm_("T", "N", x1.size(), x0.size()*a2.size()*c1.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), c1.size(), a2.size(), x0.size());
  out()->add_block(odata, x0, a2, c1, x1);
}

void Task392::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I371
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x3, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x3, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3, x2, x1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, x3, x2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), x3.size(), x2.size());
      // tensor label: I418
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), x1.size());
      dgemm_("T", "N", x3.size()*x2.size(), x0.size()*x1.size(), c3.size()*a2.size(),
             1.0, i0data_sorted, c3.size()*a2.size(), i1data_sorted, c3.size()*a2.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x0.size(), x1.size());
  out()->add_block(odata, x0, x3, x2, x1);
}

void Task393::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index x1 = b(3);
  // tensor label: I418
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, c3, a2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, c3, a2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3, a2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3, a2, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), x1.size());
    // tensor label: l2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x0)]);
    sort_indices<0,1,2,3,0,1,-1,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
    dgemm_("T", "N", x1.size(), x0.size()*c3.size()*a2.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), a2.size(), c3.size(), x0.size());
  out()->add_block(odata, x0, c3, a2, x1);
}

void Task394::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I371
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x3, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x3, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3, x2, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x3, x2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x3.size(), x2.size());
      // tensor label: I422
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a2, c1, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a2, c1, x1)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a2.size(), c1.size(), x1.size());
      dgemm_("T", "N", x3.size()*x2.size(), x0.size()*x1.size(), a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x0.size(), x1.size());
  out()->add_block(odata, x0, x3, x2, x1);
}

void Task395::Task_local::compute() {
  const Index x0 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x1 = b(3);
  // tensor label: I422
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, a2, c1, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, a2, c1, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a2, c1, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a2, c1, x1), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), x1.size());
    // tensor label: l2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x0)]);
    sort_indices<2,0,1,3,0,1,2,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
    dgemm_("T", "N", x1.size(), x0.size()*a2.size()*c1.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), c1.size(), a2.size(), x0.size());
  out()->add_block(odata, x0, a2, c1, x1);
}

void Task396::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I371
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x3, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x3, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3, x2, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: l2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c1, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, c1, x1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), c1.size(), x1.size());
    // tensor label: I714
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), x0.size());
    dgemm_("T", "N", x1.size()*x2.size()*x3.size(), x0.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x1.size()*x2.size()*x3.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size());
  out()->add_block(odata, x0, x3, x2, x1);
}

void Task397::Task_local::compute() {
  const Index c1 = b(0);
  const Index x0 = b(1);
  // tensor label: I714
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x0), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(a3, c2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, c2)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), c2.size());
      // tensor label: I715
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a3, c1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a3, c1, x0)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), a3.size(), c1.size(), x0.size());
      dgemm_("T", "N", 1, c1.size()*x0.size(), c2.size()*a3.size(),
             1.0, i0data_sorted, c2.size()*a3.size(), i1data_sorted, c2.size()*a3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x0.size());
  out()->add_block(odata, c1, x0);
}

void Task398::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: I715
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, a3, c1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c2, a3, c1, x0), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c1, x0);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, c2.size(), a3.size(), c1.size(), x0.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a3, c2, x0);
    sort_indices<2,1,0,3,1,1,-1,1>(i1data, odata, c1.size(), a3.size(), c2.size(), x0.size());
  }
  out()->add_block(odata, c2, a3, c1, x0);
}

void Task399::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I371
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x3, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x3, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3, x2, x1), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: l2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, x2)]);
      sort_indices<2,1,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), x2.size());
      // tensor label: I804
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x0.size(), x1.size());
      dgemm_("T", "N", x2.size()*x3.size(), x0.size()*x1.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, x2.size()*x3.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x0.size(), x1.size());
  out()->add_block(odata, x0, x3, x2, x1);
}

#endif
