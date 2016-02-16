//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_tasks14.cc
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

#include <src/smith/mrci/MRCI_tasks14.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task650::Task_local::compute() {
  const Index x0 = b(0);
  const Index c5 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I1135
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, c5, c3, a2)]);
  std::fill_n(odata.get(), out()->get_size(x0, c5, c3, a2), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, c5, c3, a2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x0.size(), c5.size(), c3.size(), a2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x0, a2, c3, c5);
    sort_indices<0,3,2,1,1,1,-2,1>(i1data, odata, x0.size(), a2.size(), c3.size(), c5.size());
  }
  out()->add_block(odata, x0, c5, c3, a2);
}

void Task651::Task_local::compute() {
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
      std::unique_ptr<double[]> i0data = in(0)->get_block(c5, a2, c1, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c5, a2, c1, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, c5.size(), a2.size(), c1.size(), x1.size());
      // tensor label: I1137
      std::unique_ptr<double[]> i1data = in(1)->get_block(c5, c3, a4, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c5, c3, a4, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, c5.size(), c3.size(), a4.size(), x1.size());
      dgemm_("T", "N", a2.size()*c1.size(), c3.size()*a4.size(), c5.size()*x1.size(),
             1.0, i0data_sorted, c5.size()*x1.size(), i1data_sorted, c5.size()*x1.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), c3.size(), a4.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task652::Task_local::compute() {
  const Index c5 = b(0);
  const Index c3 = b(1);
  const Index a4 = b(2);
  const Index x1 = b(3);
  // tensor label: I1137
  std::unique_ptr<double[]> odata(new double[out()->get_size(c5, c3, a4, x1)]);
  std::fill_n(odata.get(), out()->get_size(c5, c3, a4, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c5, c3, a4, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c5, c3, a4, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I1138
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c5, c3, a4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c5, c3, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c5.size(), c3.size(), a4.size());
    dgemm_("T", "N", x1.size(), c5.size()*c3.size()*a4.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), c5.size(), c3.size(), a4.size());
  out()->add_block(odata, c5, c3, a4, x1);
}

void Task653::Task_local::compute() {
  const Index x0 = b(0);
  const Index c5 = b(1);
  const Index c3 = b(2);
  const Index a4 = b(3);
  // tensor label: I1138
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, c5, c3, a4)]);
  std::fill_n(odata.get(), out()->get_size(x0, c5, c3, a4), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, c5, c3, a4);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x0.size(), c5.size(), c3.size(), a4.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x0, a4, c3, c5);
    sort_indices<0,3,2,1,1,1,1,1>(i1data, odata, x0.size(), a4.size(), c3.size(), c5.size());
  }
  out()->add_block(odata, x0, c5, c3, a4);
}

void Task654::Task_local::compute() {
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
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c5, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c5, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c5.size(), x1.size());
      // tensor label: I1146
      std::unique_ptr<double[]> i1data = in(1)->get_block(c5, c3, a2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c5, c3, a2, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, c5.size(), c3.size(), a2.size(), x1.size());
      dgemm_("T", "N", c1.size()*a4.size(), c3.size()*a2.size(), c5.size()*x1.size(),
             1.0, i0data_sorted, c5.size()*x1.size(), i1data_sorted, c5.size()*x1.size(),
             1.0, odata_sorted, c1.size()*a4.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task655::Task_local::compute() {
  const Index c5 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1146
  std::unique_ptr<double[]> odata(new double[out()->get_size(c5, c3, a2, x1)]);
  std::fill_n(odata.get(), out()->get_size(c5, c3, a2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c5, c3, a2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c5, c3, a2, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I1147
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

void Task656::Task_local::compute() {
  const Index x0 = b(0);
  const Index c5 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I1147
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, c5, c3, a2)]);
  std::fill_n(odata.get(), out()->get_size(x0, c5, c3, a2), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, c5, c3, a2);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x0.size(), c5.size(), c3.size(), a2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x0, a2, c3, c5);
    sort_indices<0,3,2,1,1,1,1,1>(i1data, odata, x0.size(), a2.size(), c3.size(), c5.size());
  }
  out()->add_block(odata, x0, c5, c3, a2);
}

void Task657::Task_local::compute() {
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
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c5, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c5, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c5.size(), x1.size());
      // tensor label: I1149
      std::unique_ptr<double[]> i1data = in(1)->get_block(c5, c3, a4, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c5, c3, a4, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, c5.size(), c3.size(), a4.size(), x1.size());
      dgemm_("T", "N", c1.size()*a2.size(), c3.size()*a4.size(), c5.size()*x1.size(),
             1.0, i0data_sorted, c5.size()*x1.size(), i1data_sorted, c5.size()*x1.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task658::Task_local::compute() {
  const Index c5 = b(0);
  const Index c3 = b(1);
  const Index a4 = b(2);
  const Index x1 = b(3);
  // tensor label: I1149
  std::unique_ptr<double[]> odata(new double[out()->get_size(c5, c3, a4, x1)]);
  std::fill_n(odata.get(), out()->get_size(c5, c3, a4, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c5, c3, a4, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c5, c3, a4, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I1150
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c5, c3, a4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c5, c3, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c5.size(), c3.size(), a4.size());
    dgemm_("T", "N", x1.size(), c5.size()*c3.size()*a4.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), c5.size(), c3.size(), a4.size());
  out()->add_block(odata, c5, c3, a4, x1);
}

void Task659::Task_local::compute() {
  const Index x0 = b(0);
  const Index c5 = b(1);
  const Index c3 = b(2);
  const Index a4 = b(3);
  // tensor label: I1150
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, c5, c3, a4)]);
  std::fill_n(odata.get(), out()->get_size(x0, c5, c3, a4), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, c5, c3, a4);
    sort_indices<0,1,2,3,1,1,4,1>(i0data, odata, x0.size(), c5.size(), c3.size(), a4.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x0, a4, c3, c5);
    sort_indices<0,3,2,1,1,1,-2,1>(i1data, odata, x0.size(), a4.size(), c3.size(), c5.size());
  }
  out()->add_block(odata, x0, c5, c3, a4);
}

void Task660::Task_local::compute() {
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
    for (auto& a5 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a4, a5, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a4, a5, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a4.size(), a5.size(), a2.size());
      // tensor label: I1158
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a5, c3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a5, c3, x0)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), a5.size(), c3.size(), x0.size());
      dgemm_("T", "N", a4.size()*a2.size(), c1.size()*c3.size(), a5.size()*x0.size(),
             1.0, i0data_sorted, a5.size()*x0.size(), i1data_sorted, a5.size()*x0.size(),
             1.0, odata_sorted, a4.size()*a2.size());
    }
  }
  sort_indices<1,2,0,3,1,1,1,1>(odata_sorted, odata, a4.size(), a2.size(), c1.size(), c3.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task661::Task_local::compute() {
  const Index c1 = b(0);
  const Index a5 = b(1);
  const Index c3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1158
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a5, c3, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, a5, c3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a5, c3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a5, c3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a5, c3, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a5, c3, x1)]);
    sort_indices<3,0,1,2,0,1,2,1>(i1data, i1data_sorted, c1.size(), a5.size(), c3.size(), x1.size());
    dgemm_("T", "N", x0.size(), c1.size()*a5.size()*c3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x0.size(), c1.size(), a5.size(), c3.size());
  out()->add_block(odata, c1, a5, c3, x0);
}

void Task662::Task_local::compute() {
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
    for (auto& a5 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a2, a5, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a2, a5, a4)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a2.size(), a5.size(), a4.size());
      // tensor label: I1161
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a5, c3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a5, c3, x0)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), a5.size(), c3.size(), x0.size());
      dgemm_("T", "N", a2.size()*a4.size(), c1.size()*c3.size(), a5.size()*x0.size(),
             1.0, i0data_sorted, a5.size()*x0.size(), i1data_sorted, a5.size()*x0.size(),
             1.0, odata_sorted, a2.size()*a4.size());
    }
  }
  sort_indices<0,2,1,3,1,1,1,1>(odata_sorted, odata, a2.size(), a4.size(), c1.size(), c3.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task663::Task_local::compute() {
  const Index c1 = b(0);
  const Index a5 = b(1);
  const Index c3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1161
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a5, c3, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, a5, c3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a5, c3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a5, c3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a5, c3, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a5, c3, x1)]);
    sort_indices<3,0,1,2,0,1,-1,1>(i1data, i1data_sorted, c1.size(), a5.size(), c3.size(), x1.size());
    dgemm_("T", "N", x0.size(), c1.size()*a5.size()*c3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x0.size(), c1.size(), a5.size(), c3.size());
  out()->add_block(odata, c1, a5, c3, x0);
}

void Task664::Task_local::compute() {
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
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a4, c1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c1.size(), x2.size());
      // tensor label: I1164
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), x3.size(), x2.size());
      dgemm_("T", "N", a4.size()*c1.size(), c3.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a4.size()*c1.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), c3.size(), a2.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task665::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1164
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a2, x3, x2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma29
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I1165
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), x1.size(), x0.size());
      dgemm_("T", "N", x3.size()*x2.size(), c3.size()*a2.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c3.size(), a2.size());
  out()->add_block(odata, c3, a2, x3, x2);
}

void Task666::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1165
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, x1, x0), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, x1, x0);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, c3.size(), a2.size(), x1.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x1, x0, c3, a2);
    sort_indices<2,3,0,1,1,1,1,2>(i1data, odata, x1.size(), x0.size(), c3.size(), a2.size());
  }
  out()->add_block(odata, c3, a2, x1, x0);
}

void Task667::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1164
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a2, x3, x2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma18
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x0, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x0.size(), x2.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x1, x0, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x1, x0, a2)]);
      sort_indices<1,2,0,3,0,1,1,2>(i1data, i1data_sorted, c3.size(), x1.size(), x0.size(), a2.size());
      dgemm_("T", "N", x3.size()*x2.size(), c3.size()*a2.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c3.size(), a2.size());
  out()->add_block(odata, c3, a2, x3, x2);
}

void Task668::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1164
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a2, x3, x2), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma27
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x1, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a2, c3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a2, c3, x0)]);
      sort_indices<3,0,1,2,0,1,-1,2>(i1data, i1data_sorted, x1.size(), a2.size(), c3.size(), x0.size());
      dgemm_("T", "N", x3.size()*x2.size(), a2.size()*c3.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<3,2,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), a2.size(), c3.size());
  out()->add_block(odata, c3, a2, x3, x2);
}

void Task669::Task_local::compute() {
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
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c1.size(), x2.size());
      // tensor label: I1167
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a4, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a4, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), a4.size(), x3.size(), x2.size());
      dgemm_("T", "N", a2.size()*c1.size(), c3.size()*a4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), c3.size(), a4.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task670::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1167
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a4, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c3, a4, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a4, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a4, x3, x2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma29
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I1168
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a4, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a4, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), a4.size(), x1.size(), x0.size());
      dgemm_("T", "N", x3.size()*x2.size(), c3.size()*a4.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c3.size(), a4.size());
  out()->add_block(odata, c3, a4, x3, x2);
}

void Task671::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1168
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a4, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, a4, x1, x0), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c3.size(), a4.size(), x1.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x1, a4, c3, x0);
    sort_indices<2,1,0,3,1,1,1,2>(i1data, odata, x1.size(), a4.size(), c3.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i2data = in(0)->get_block(x1, x0, c3, a4);
    sort_indices<2,3,0,1,1,1,-1,1>(i2data, odata, x1.size(), x0.size(), c3.size(), a4.size());
  }
  out()->add_block(odata, c3, a4, x1, x0);
}

void Task672::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1167
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a4, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c3, a4, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a4, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a4, x3, x2), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma10
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x1, x0, a4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x1, x0, a4)]);
      sort_indices<2,1,0,3,0,1,-1,2>(i1data, i1data_sorted, c3.size(), x1.size(), x0.size(), a4.size());
      dgemm_("T", "N", x3.size()*x2.size(), c3.size()*a4.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c3.size(), a4.size());
  out()->add_block(odata, c3, a4, x3, x2);
}

void Task673::Task_local::compute() {
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
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, c5);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, c5)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), c5.size());
    // tensor label: I1188
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c5);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c5)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a4.size(), c5.size());
    dgemm_("T", "N", c1.size()*a2.size()*c3.size(), a4.size(), c5.size(),
           1.0, i0data_sorted, c5.size(), i1data_sorted, c5.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task674::Task_local::compute() {
  const Index a4 = b(0);
  const Index c5 = b(1);
  // tensor label: I1188
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, c5)]);
  std::fill_n(odata.get(), out()->get_size(a4, c5), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c5), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma32
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: I1189
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a4, c5, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a4, c5, x0)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), a4.size(), c5.size(), x0.size());
      dgemm_("T", "N", 1, a4.size()*c5.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a4.size(), c5.size());
  out()->add_block(odata, a4, c5);
}

void Task675::Task_local::compute() {
  const Index x1 = b(0);
  const Index a4 = b(1);
  const Index c5 = b(2);
  const Index x0 = b(3);
  // tensor label: I1189
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, a4, c5, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, a4, c5, x0), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c5, x0);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x1.size(), a4.size(), c5.size(), x0.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c5, a4, x1, x0);
    sort_indices<2,1,0,3,1,1,-4,1>(i1data, odata, c5.size(), a4.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x1, a4, c5, x0);
}

void Task676::Task_local::compute() {
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
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, c5);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, c5)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), c5.size());
    // tensor label: I1191
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c5);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c5)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a2.size(), c5.size());
    dgemm_("T", "N", c1.size()*a4.size()*c3.size(), a2.size(), c5.size(),
           1.0, i0data_sorted, c5.size(), i1data_sorted, c5.size(),
           1.0, odata_sorted, c1.size()*a4.size()*c3.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task677::Task_local::compute() {
  const Index a2 = b(0);
  const Index c5 = b(1);
  // tensor label: I1191
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c5)]);
  std::fill_n(odata.get(), out()->get_size(a2, c5), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c5), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma32
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: I1192
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a2, c5, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a2, c5, x0)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), a2.size(), c5.size(), x0.size());
      dgemm_("T", "N", 1, a2.size()*c5.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a2.size(), c5.size());
  out()->add_block(odata, a2, c5);
}

void Task678::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  const Index c5 = b(2);
  const Index x0 = b(3);
  // tensor label: I1192
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, a2, c5, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, a2, c5, x0), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c5, x0);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x1.size(), a2.size(), c5.size(), x0.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c5, a2, x1, x0);
    sort_indices<2,1,0,3,1,1,2,1>(i1data, odata, c5.size(), a2.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x1, a2, c5, x0);
}

void Task679::Task_local::compute() {
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
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, a5, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4, a5, a2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), a5.size(), a2.size());
    // tensor label: I1194
    std::unique_ptr<double[]> i1data = in(1)->get_block(a5, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, c1)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a5.size(), c1.size());
    dgemm_("T", "N", c3.size()*a4.size()*a2.size(), c1.size(), a5.size(),
           1.0, i0data_sorted, a5.size(), i1data_sorted, a5.size(),
           1.0, odata_sorted, c3.size()*a4.size()*a2.size());
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, c3.size(), a4.size(), a2.size(), c1.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task680::Task_local::compute() {
  const Index a5 = b(0);
  const Index c1 = b(1);
  // tensor label: I1194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a5, c1)]);
  std::fill_n(odata.get(), out()->get_size(a5, c1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a5, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a5, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma32
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: I1195
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a5, c1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a5, c1, x0)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), a5.size(), c1.size(), x0.size());
      dgemm_("T", "N", 1, a5.size()*c1.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a5.size(), c1.size());
  out()->add_block(odata, a5, c1);
}

void Task681::Task_local::compute() {
  const Index x1 = b(0);
  const Index a5 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1195
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, a5, c1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, a5, c1, x0), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a5, c1, x0);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x1.size(), a5.size(), c1.size(), x0.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a5, x1, x0);
    sort_indices<2,1,0,3,1,1,4,1>(i1data, odata, c1.size(), a5.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x1, a5, c1, x0);
}

void Task682::Task_local::compute() {
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
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, a5, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, a5, a4)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), a5.size(), a4.size());
    // tensor label: I1197
    std::unique_ptr<double[]> i1data = in(1)->get_block(a5, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, c1)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a5.size(), c1.size());
    dgemm_("T", "N", c3.size()*a2.size()*a4.size(), c1.size(), a5.size(),
           1.0, i0data_sorted, a5.size(), i1data_sorted, a5.size(),
           1.0, odata_sorted, c3.size()*a2.size()*a4.size());
  }
  sort_indices<1,3,2,0,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), a4.size(), c1.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task683::Task_local::compute() {
  const Index a5 = b(0);
  const Index c1 = b(1);
  // tensor label: I1197
  std::unique_ptr<double[]> odata(new double[out()->get_size(a5, c1)]);
  std::fill_n(odata.get(), out()->get_size(a5, c1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a5, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a5, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma32
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: I1198
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a5, c1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a5, c1, x0)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), a5.size(), c1.size(), x0.size());
      dgemm_("T", "N", 1, a5.size()*c1.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a5.size(), c1.size());
  out()->add_block(odata, a5, c1);
}

void Task684::Task_local::compute() {
  const Index x1 = b(0);
  const Index a5 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1198
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, a5, c1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, a5, c1, x0), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a5, c1, x0);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x1.size(), a5.size(), c1.size(), x0.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a5, x1, x0);
    sort_indices<2,1,0,3,1,1,-2,1>(i1data, odata, c1.size(), a5.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x1, a5, c1, x0);
}

void Task685::Task_local::compute() {
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
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), x3.size(), x2.size());
      // tensor label: I1200
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), x3.size(), x2.size());
      dgemm_("T", "N", c1.size()*a4.size(), c3.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c1.size()*a4.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task686::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1200
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a2, x3, x2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma29
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I1201
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), x1.size(), x0.size());
      dgemm_("T", "N", x3.size()*x2.size(), c3.size()*a2.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c3.size(), a2.size());
  out()->add_block(odata, c3, a2, x3, x2);
}

void Task687::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1201
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, x1, x0), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c3.size(), a2.size(), x1.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x1, a2, c3, x0);
    sort_indices<2,1,0,3,1,1,1,2>(i1data, odata, x1.size(), a2.size(), c3.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i2data = in(0)->get_block(x1, x0, c3, a2);
    sort_indices<2,3,0,1,1,1,-1,1>(i2data, odata, x1.size(), x0.size(), c3.size(), a2.size());
  }
  out()->add_block(odata, c3, a2, x1, x0);
}

void Task688::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1200
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a2, x3, x2), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma10
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x1, x0, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x1, x0, a2)]);
      sort_indices<2,1,0,3,0,1,-1,2>(i1data, i1data_sorted, c3.size(), x1.size(), x0.size(), a2.size());
      dgemm_("T", "N", x3.size()*x2.size(), c3.size()*a2.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c3.size(), a2.size());
  out()->add_block(odata, c3, a2, x3, x2);
}

void Task689::Task_local::compute() {
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
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x3.size(), x2.size());
      // tensor label: I1203
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a4, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a4, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), a4.size(), x3.size(), x2.size());
      dgemm_("T", "N", c1.size()*a2.size(), c3.size()*a4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task690::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1203
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a4, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c3, a4, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a4, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a4, x3, x2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma29
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I1204
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a4, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a4, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), a4.size(), x1.size(), x0.size());
      dgemm_("T", "N", x3.size()*x2.size(), c3.size()*a4.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c3.size(), a4.size());
  out()->add_block(odata, c3, a4, x3, x2);
}

void Task691::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1204
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a4, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, a4, x1, x0), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, x1, x0);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, c3.size(), a4.size(), x1.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x1, a4, c3, x0);
    sort_indices<2,1,0,3,1,1,-1,1>(i1data, odata, x1.size(), a4.size(), c3.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i2data = in(0)->get_block(x1, x0, c3, a4);
    sort_indices<2,3,0,1,1,1,2,1>(i2data, odata, x1.size(), x0.size(), c3.size(), a4.size());
  }
  out()->add_block(odata, c3, a4, x1, x0);
}

void Task692::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1203
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a4, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c3, a4, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a4, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a4, x3, x2), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma10
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x1, x0, a4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x1, x0, a4)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), x1.size(), x0.size(), a4.size());
      dgemm_("T", "N", x3.size()*x2.size(), c3.size()*a4.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c3.size(), a4.size());
  out()->add_block(odata, c3, a4, x3, x2);
}

void Task693::Task_local::compute() {
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
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x0, c3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x0, c3, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x0.size(), c3.size(), a2.size());
    // tensor label: I1236
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x0)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a4.size(), x0.size());
    dgemm_("T", "N", c1.size()*c3.size()*a2.size(), a4.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, c1.size()*c3.size()*a2.size());
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), a2.size(), a4.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task694::Task_local::compute() {
  const Index a4 = b(0);
  const Index x0 = b(1);
  // tensor label: I1236
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma51
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a4, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a4, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), a4.size(), x2.size(), x1.size());
        dgemm_("T", "N", x0.size(), a4.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a4.size());
  out()->add_block(odata, a4, x0);
}

void Task695::Task_local::compute() {
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
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x0, c3, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x0, c3, a4)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x0.size(), c3.size(), a4.size());
    // tensor label: I1239
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x0)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a2.size(), x0.size());
    dgemm_("T", "N", c1.size()*c3.size()*a4.size(), a2.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, c1.size()*c3.size()*a4.size());
  }
  sort_indices<3,0,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), a4.size(), a2.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task696::Task_local::compute() {
  const Index a2 = b(0);
  const Index x0 = b(1);
  // tensor label: I1239
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, x0)]);
  std::fill_n(odata.get(), out()->get_size(a2, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma51
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a2, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a2, x2, x1)]);
        sort_indices<0,2,3,1,0,1,-2,1>(i1data, i1data_sorted, x3.size(), a2.size(), x2.size(), x1.size());
        dgemm_("T", "N", x0.size(), a2.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a2.size());
  out()->add_block(odata, a2, x0);
}

void Task697::Task_local::compute() {
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
    for (auto& c6 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a5, c6, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a5, c6, a4)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a5.size(), c6.size(), a4.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, c6, a5, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, c6, a5, a2)]);
      sort_indices<2,1,0,3,0,1,-8,1>(i1data, i1data_sorted, c3.size(), c6.size(), a5.size(), a2.size());
      dgemm_("T", "N", c1.size()*a4.size(), c3.size()*a2.size(), c6.size()*a5.size(),
             1.0, i0data_sorted, c6.size()*a5.size(), i1data_sorted, c6.size()*a5.size(),
             1.0, odata_sorted, c1.size()*a4.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task698::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& c6 : *range_[0]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c6, a5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c6, a5)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c6.size(), a5.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, c6, a5, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, c6, a5, a2)]);
      sort_indices<1,2,0,3,0,1,4,1>(i1data, i1data_sorted, c3.size(), c6.size(), a5.size(), a2.size());
      dgemm_("T", "N", c1.size()*a4.size(), c3.size()*a2.size(), c6.size()*a5.size(),
             1.0, i0data_sorted, c6.size()*a5.size(), i1data_sorted, c6.size()*a5.size(),
             1.0, odata_sorted, c1.size()*a4.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

void Task699::Task_local::compute() {
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
    for (auto& c6 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a5, c6, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a5, c6, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a5.size(), c6.size(), a2.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, c6, a5, a4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, c6, a5, a4)]);
      sort_indices<2,1,0,3,0,1,4,1>(i1data, i1data_sorted, c3.size(), c6.size(), a5.size(), a4.size());
      dgemm_("T", "N", c1.size()*a2.size(), c3.size()*a4.size(), c6.size()*a5.size(),
             1.0, i0data_sorted, c6.size()*a5.size(), i1data_sorted, c6.size()*a5.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->add_block(odata, a2, c1, a4, c3);
}

#endif
