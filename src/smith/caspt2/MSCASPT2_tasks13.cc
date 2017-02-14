//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_tasks13.cc
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

#include <src/smith/caspt2/MSCASPT2_tasks13.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MSCASPT2;

void Task600::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I868
  std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
  // tensor label (calculated on-the-fly): Gamma148
  {
    std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x1, x0)]);
    std::fill_n(o1data.get(), out(1)->get_size(x1, x0), 0.0);
    sort_indices<0,1,1,1,1,1>(i0data_sorted, o1data, x1.size(), x0.size());
    out(1)->add_block(o1data, x1, x0);
  }
}

void Task601::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I868
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c2, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, c2, a1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c2.size(), a1.size());
        // tensor label: l2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a1, c2, a3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a1, c2, a3)]);
        sort_indices<3,2,1,0,0,1,-1,1>(i1data, i1data_sorted, x1.size(), a1.size(), c2.size(), a3.size());
        dgemm_("T", "N", x0.size(), x1.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->add_block(odata, x1, x0);
}

void Task602::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I868
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
        // tensor label: l2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a1, c2, a3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a1, c2, a3)]);
        sort_indices<1,2,3,0,0,1,2,1>(i1data, i1data_sorted, x1.size(), a1.size(), c2.size(), a3.size());
        dgemm_("T", "N", x0.size(), x1.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->add_block(odata, x1, x0);
}

void Task603::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I868
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: h1
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size());
      // tensor label: l2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a1, c2, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a1, c2, x0)]);
      sort_indices<2,1,0,3,0,1,-1,1>(i1data, i1data_sorted, x1.size(), a1.size(), c2.size(), x0.size());
      dgemm_("T", "N", 1, x0.size()*x1.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->add_block(odata, x1, x0);
}

void Task604::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I868
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: h1
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size());
      // tensor label: l2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,2,1>(i1data, i1data_sorted, c1.size(), a2.size(), x1.size(), x0.size());
      dgemm_("T", "N", 1, x0.size()*x1.size(), a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->add_block(odata, x1, x0);
}

void Task605::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I874
  std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x0, x1);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x3, x0, x1)]);
  sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x0.size(), x1.size());
  // tensor label (calculated on-the-fly): Gamma170
  {
    std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x3, x0, x2, x1)]);
    std::fill_n(o2data.get(), out(2)->get_size(x3, x0, x2, x1), 0.0);
    sort_indices<0,1,2,3,1,1,1,1>(i0data_sorted, o2data, x3.size(), x0.size(), x2.size(), x1.size());
    out(2)->add_block(o2data, x3, x0, x2, x1);
  }
}

void Task606::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I874
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, x3, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x2, x3, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x3, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, x0, x1), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());
      // tensor label: l2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, x2, a2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), x2.size(), a2.size());
      dgemm_("T", "N", x0.size()*x1.size(), x2.size()*x3.size(), a2.size()*a1.size(),
             1.0, i0data_sorted, a2.size()*a1.size(), i1data_sorted, a2.size()*a1.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<3,2,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x3.size(), x2.size());
  out()->add_block(odata, x2, x3, x0, x1);
}

void Task607::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I874
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, x3, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x2, x3, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x3, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, x0, x1), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size());
    // tensor label: l2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, x2, x1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), x2.size(), x1.size());
    dgemm_("T", "N", x0.size(), x1.size()*x2.size()*x3.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<2,1,0,3,1,1,1,1>(odata_sorted, odata, x0.size(), x3.size(), x2.size(), x1.size());
  out()->add_block(odata, x2, x3, x0, x1);
}

#endif
