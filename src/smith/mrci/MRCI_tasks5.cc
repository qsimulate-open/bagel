//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_tasks5.cc
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

#include <src/smith/mrci/MRCI_tasks5.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task200::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I37
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x1.size());
    dgemm_("T", "N", x0.size(), c1.size()*a2.size()*c3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,3,0,2,1,1,1,1>(odata_sorted, odata, x0.size(), c1.size(), a2.size(), c3.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task201::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a2, c3, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c3.size(), x1.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c4)]);
    sort_indices<1,0,0,1,-2,1>(i1data, i1data_sorted, c1.size(), c4.size());
    dgemm_("T", "N", a2.size()*c3.size()*x1.size(), c1.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, a2.size()*c3.size()*x1.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a2.size(), c3.size(), x1.size(), c1.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task202::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c4, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c4, x1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c4.size(), x1.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c4)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), c4.size());
    dgemm_("T", "N", c3.size()*a2.size()*x1.size(), c1.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c3.size()*a2.size()*x1.size());
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), x1.size(), c1.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task203::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c1, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a2, c1, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c1.size(), x1.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, c4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, c4)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), c4.size());
    dgemm_("T", "N", a2.size()*c1.size()*x1.size(), c3.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, a2.size()*c1.size()*x1.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x1.size(), c3.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task204::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, c1, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4, c1, x1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), c1.size(), x1.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a2)]);
    sort_indices<0,1,0,1,-1,1>(i1data, i1data_sorted, a4.size(), a2.size());
    dgemm_("T", "N", c3.size()*c1.size()*x1.size(), a2.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, c3.size()*c1.size()*x1.size());
  }
  sort_indices<1,3,0,2,1,1,1,1>(odata_sorted, odata, c3.size(), c1.size(), x1.size(), a2.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task205::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c4, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c4, x1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c4.size(), x1.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, c4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, c4)]);
    sort_indices<1,0,0,1,-2,1>(i1data, i1data_sorted, c3.size(), c4.size());
    dgemm_("T", "N", c1.size()*a2.size()*x1.size(), c3.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c1.size()*a2.size()*x1.size());
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), c3.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task206::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, x1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x1.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a2)]);
    sort_indices<0,1,0,1,2,1>(i1data, i1data_sorted, a4.size(), a2.size());
    dgemm_("T", "N", c1.size()*c3.size()*x1.size(), a2.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, c1.size()*c3.size()*x1.size());
  }
  sort_indices<0,3,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), x1.size(), a2.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task207::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c5 : *range_[0]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c5, a2, c4, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c5, a2, c4, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c5.size(), a2.size(), c4.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c5, c3, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c5, c3, c4)]);
      sort_indices<1,3,0,2,0,1,2,1>(i1data, i1data_sorted, c1.size(), c5.size(), c3.size(), c4.size());
      dgemm_("T", "N", a2.size()*x1.size(), c1.size()*c3.size(), c5.size()*c4.size(),
             1.0, i0data_sorted, c5.size()*c4.size(), i1data_sorted, c5.size()*c4.size(),
             1.0, odata_sorted, a2.size()*x1.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, a2.size(), x1.size(), c1.size(), c3.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task208::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& c5 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c5, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a2, c5, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c5.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c5, c3, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c5, c3, c4)]);
      sort_indices<3,1,0,2,0,1,-1,1>(i1data, i1data_sorted, c1.size(), c5.size(), c3.size(), c4.size());
      dgemm_("T", "N", a2.size()*x1.size(), c1.size()*c3.size(), c5.size()*c4.size(),
             1.0, i0data_sorted, c5.size()*c4.size(), i1data_sorted, c5.size()*c4.size(),
             1.0, odata_sorted, a2.size()*x1.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, a2.size(), x1.size(), c1.size(), c3.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task209::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c5 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c5, a4, c3, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c5, a4, c3, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c5.size(), a4.size(), c3.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c5, a4, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c5, a4, a2)]);
      sort_indices<1,2,0,3,0,1,-2,1>(i1data, i1data_sorted, c1.size(), c5.size(), a4.size(), a2.size());
      dgemm_("T", "N", c3.size()*x1.size(), c1.size()*a2.size(), c5.size()*a4.size(),
             1.0, i0data_sorted, c5.size()*a4.size(), i1data_sorted, c5.size()*a4.size(),
             1.0, odata_sorted, c3.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c3.size(), x1.size(), c1.size(), a2.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task210::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a5, c3, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a5, c3, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a5.size(), c3.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, a5, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, a5, c4)]);
      sort_indices<3,2,0,1,0,1,4,1>(i1data, i1data_sorted, c1.size(), a2.size(), a5.size(), c4.size());
      dgemm_("T", "N", c3.size()*x1.size(), c1.size()*a2.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, c3.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c3.size(), x1.size(), c1.size(), a2.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task211::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c5 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, c5, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4, c5, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), c5.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c5, a4, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c5, a4, a2)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), c5.size(), a4.size(), a2.size());
      dgemm_("T", "N", c3.size()*x1.size(), c1.size()*a2.size(), c5.size()*a4.size(),
             1.0, i0data_sorted, c5.size()*a4.size(), i1data_sorted, c5.size()*a4.size(),
             1.0, odata_sorted, c3.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c3.size(), x1.size(), c1.size(), a2.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task212::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a5 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a5, c4, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a5, c4, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a5.size(), c4.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, a5, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, a5, c4)]);
      sort_indices<2,3,0,1,0,1,-2,1>(i1data, i1data_sorted, c1.size(), a2.size(), a5.size(), c4.size());
      dgemm_("T", "N", c3.size()*x1.size(), c1.size()*a2.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, c3.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c3.size(), x1.size(), c1.size(), a2.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task213::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c5 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c5, a4, c1, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c5, a4, c1, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c5.size(), a4.size(), c1.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, c5, a4, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, c5, a4, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), c5.size(), a4.size(), a2.size());
      dgemm_("T", "N", c1.size()*x1.size(), c3.size()*a2.size(), c5.size()*a4.size(),
             1.0, i0data_sorted, c5.size()*a4.size(), i1data_sorted, c5.size()*a4.size(),
             1.0, odata_sorted, c1.size()*x1.size());
    }
  }
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), x1.size(), c3.size(), a2.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task214::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a5, c1, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a5, c1, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a5.size(), c1.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, a5, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, a5, c4)]);
      sort_indices<3,2,0,1,0,1,-2,1>(i1data, i1data_sorted, c3.size(), a2.size(), a5.size(), c4.size());
      dgemm_("T", "N", c1.size()*x1.size(), c3.size()*a2.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, c1.size()*x1.size());
    }
  }
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), x1.size(), c3.size(), a2.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task215::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c5 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c5, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c5, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c5.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, c5, a4, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, c5, a4, a2)]);
      sort_indices<2,1,0,3,0,1,-2,1>(i1data, i1data_sorted, c3.size(), c5.size(), a4.size(), a2.size());
      dgemm_("T", "N", c1.size()*x1.size(), c3.size()*a2.size(), c5.size()*a4.size(),
             1.0, i0data_sorted, c5.size()*a4.size(), i1data_sorted, c5.size()*a4.size(),
             1.0, odata_sorted, c1.size()*x1.size());
    }
  }
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), x1.size(), c3.size(), a2.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task216::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a5 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a5, c4, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a5, c4, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a5.size(), c4.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, a5, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, a5, c4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), a5.size(), c4.size());
      dgemm_("T", "N", c1.size()*x1.size(), c3.size()*a2.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, c1.size()*x1.size());
    }
  }
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), x1.size(), c3.size(), a2.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task217::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a5 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a5, c4, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a5, c4, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a5.size(), c4.size(), a2.size());
      // tensor label: I613
      std::unique_ptr<double[]> i1data = in(1)->get_block(a5, x1, c1, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, x1, c1, c4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, a5.size(), x1.size(), c1.size(), c4.size());
      dgemm_("T", "N", c3.size()*a2.size(), x1.size()*c1.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, c3.size()*a2.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), x1.size(), c1.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task218::Task_local::compute() {
  const Index a5 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index c4 = b(3);
  // tensor label: I613
  std::unique_ptr<double[]> odata(new double[out()->get_size(a5, x1, c1, c4)]);
  std::fill_n(odata.get(), out()->get_size(a5, x1, c1, c4), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, x1, c1, c4);
    sort_indices<0,1,2,3,1,1,-4,1>(i0data, odata, a5.size(), x1.size(), c1.size(), c4.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, x1, a5, c4);
    sort_indices<2,1,0,3,1,1,2,1>(i1data, odata, c1.size(), x1.size(), a5.size(), c4.size());
  }
  out()->add_block(odata, a5, x1, c1, c4);
}

void Task219::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c4, a5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c4, a5)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c4.size(), a5.size());
      // tensor label: I616
      std::unique_ptr<double[]> i1data = in(1)->get_block(a5, x1, c1, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, x1, c1, c4)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a5.size(), x1.size(), c1.size(), c4.size());
      dgemm_("T", "N", c3.size()*a2.size(), x1.size()*c1.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, c3.size()*a2.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), x1.size(), c1.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task220::Task_local::compute() {
  const Index a5 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index c4 = b(3);
  // tensor label: I616
  std::unique_ptr<double[]> odata(new double[out()->get_size(a5, x1, c1, c4)]);
  std::fill_n(odata.get(), out()->get_size(a5, x1, c1, c4), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, x1, c1, c4);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, a5.size(), x1.size(), c1.size(), c4.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, x1, a5, c4);
    sort_indices<2,1,0,3,1,1,-4,1>(i1data, odata, c1.size(), x1.size(), a5.size(), c4.size());
  }
  out()->add_block(odata, a5, x1, c1, c4);
}

void Task221::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a5 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a5, c4, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a5, c4, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a5.size(), c4.size(), a2.size());
      // tensor label: I625
      std::unique_ptr<double[]> i1data = in(1)->get_block(a5, x1, c3, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, x1, c3, c4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, a5.size(), x1.size(), c3.size(), c4.size());
      dgemm_("T", "N", c1.size()*a2.size(), x1.size()*c3.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), c3.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task222::Task_local::compute() {
  const Index a5 = b(0);
  const Index x1 = b(1);
  const Index c3 = b(2);
  const Index c4 = b(3);
  // tensor label: I625
  std::unique_ptr<double[]> odata(new double[out()->get_size(a5, x1, c3, c4)]);
  std::fill_n(odata.get(), out()->get_size(a5, x1, c3, c4), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, x1, c3, c4);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, a5.size(), x1.size(), c3.size(), c4.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c3, x1, a5, c4);
    sort_indices<2,1,0,3,1,1,-4,1>(i1data, odata, c3.size(), x1.size(), a5.size(), c4.size());
  }
  out()->add_block(odata, a5, x1, c3, c4);
}

void Task223::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c4, a5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c4, a5)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c4.size(), a5.size());
      // tensor label: I628
      std::unique_ptr<double[]> i1data = in(1)->get_block(a5, x1, c3, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, x1, c3, c4)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a5.size(), x1.size(), c3.size(), c4.size());
      dgemm_("T", "N", c1.size()*a2.size(), x1.size()*c3.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), c3.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task224::Task_local::compute() {
  const Index a5 = b(0);
  const Index x1 = b(1);
  const Index c3 = b(2);
  const Index c4 = b(3);
  // tensor label: I628
  std::unique_ptr<double[]> odata(new double[out()->get_size(a5, x1, c3, c4)]);
  std::fill_n(odata.get(), out()->get_size(a5, x1, c3, c4), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, x1, c3, c4);
    sort_indices<0,1,2,3,1,1,-4,1>(i0data, odata, a5.size(), x1.size(), c3.size(), c4.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c3, x1, a5, c4);
    sort_indices<2,1,0,3,1,1,8,1>(i1data, odata, c3.size(), x1.size(), a5.size(), c4.size());
  }
  out()->add_block(odata, a5, x1, c3, c4);
}

void Task225::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a5 : *range_[2]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a5, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a5, c3, a4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a5.size(), c3.size(), a4.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(a5, x1, a4, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, x1, a4, a2)]);
      sort_indices<0,2,1,3,0,1,-2,1>(i1data, i1data_sorted, a5.size(), x1.size(), a4.size(), a2.size());
      dgemm_("T", "N", c1.size()*c3.size(), x1.size()*a2.size(), a5.size()*a4.size(),
             1.0, i0data_sorted, a5.size()*a4.size(), i1data_sorted, a5.size()*a4.size(),
             1.0, odata_sorted, c1.size()*c3.size());
    }
  }
  sort_indices<0,3,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), x1.size(), a2.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task226::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a5)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a5.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(a5, x1, a4, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, x1, a4, a2)]);
      sort_indices<2,0,1,3,0,1,4,1>(i1data, i1data_sorted, a5.size(), x1.size(), a4.size(), a2.size());
      dgemm_("T", "N", c1.size()*c3.size(), x1.size()*a2.size(), a5.size()*a4.size(),
             1.0, i0data_sorted, a5.size()*a4.size(), i1data_sorted, a5.size()*a4.size(),
             1.0, odata_sorted, c1.size()*c3.size());
    }
  }
  sort_indices<0,3,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), x1.size(), a2.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task227::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), x1.size());
    // tensor label: I55
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c3, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c3, x1, x0)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), c3.size(), x1.size(), x0.size());
    dgemm_("T", "N", c1.size(), a2.size()*c3.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size());
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task228::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I55
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c3, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a2, c3, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma18
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x0, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x0.size(), x2.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a2, c3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a2, c3, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a2.size(), c3.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), a2.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a2.size(), c3.size());
  out()->add_block(odata, a2, c3, x1, x0);
}

void Task229::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I55
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c3, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a2, c3, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma10
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,-1,1>(i1data, i1data_sorted, c3.size(), a2.size(), x3.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), c3.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c3.size(), a2.size());
  out()->add_block(odata, a2, c3, x1, x0);
}

void Task230::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), x1.size());
    // tensor label: I58
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c1, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c1, x0, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a2.size(), c1.size(), x0.size(), x1.size());
    dgemm_("T", "N", c3.size(), a2.size()*c1.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c3.size());
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task231::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I58
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma10
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: I59
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a2, c1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a2, c1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a2.size(), c1.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), a2.size()*c1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a2.size(), c1.size());
  out()->add_block(odata, a2, c1, x0, x1);
}

void Task232::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x2 = b(3);
  // tensor label: I59
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a2, c1, x2)]);
  std::fill_n(odata.get(), out()->get_size(x3, a2, c1, x2), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), a2.size(), c1.size(), x2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, x3, x2);
    sort_indices<2,1,0,3,1,1,2,1>(i1data, odata, c1.size(), a2.size(), x3.size(), x2.size());
  }
  out()->add_block(odata, x3, a2, c1, x2);
}

void Task233::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
    // tensor label: I67
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), x0.size());
    dgemm_("T", "N", c1.size()*c3.size()*a2.size(), x0.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, c1.size()*c3.size()*a2.size());
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), a2.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task234::Task_local::compute() {
  const Index a4 = b(0);
  const Index x0 = b(1);
  // tensor label: I67
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x1)]);
    sort_indices<1,0,0,1,-2,1>(i1data, i1data_sorted, a4.size(), x1.size());
    dgemm_("T", "N", x0.size(), a4.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a4.size());
  out()->add_block(odata, a4, x0);
}

void Task235::Task_local::compute() {
  const Index a4 = b(0);
  const Index x0 = b(1);
  // tensor label: I67
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma197
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x3, x2, x1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x3, x2, x1)]);
        sort_indices<1,2,3,0,0,1,-1,1>(i1data, i1data_sorted, a4.size(), x3.size(), x2.size(), x1.size());
        dgemm_("T", "N", x0.size(), a4.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a4.size());
  out()->add_block(odata, a4, x0);
}

void Task236::Task_local::compute() {
  const Index a4 = b(0);
  const Index x0 = b(1);
  // tensor label: I67
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma10
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, a4, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, a4, x1)]);
        sort_indices<0,1,3,2,0,1,-1,1>(i1data, i1data_sorted, x3.size(), x2.size(), a4.size(), x1.size());
        dgemm_("T", "N", x0.size(), a4.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a4.size());
  out()->add_block(odata, a4, x0);
}

void Task237::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
    // tensor label: I70
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), x0.size());
    dgemm_("T", "N", c1.size()*a2.size()*c3.size(), x0.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task238::Task_local::compute() {
  const Index a4 = b(0);
  const Index x0 = b(1);
  // tensor label: I70
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x1)]);
    sort_indices<1,0,0,1,4,1>(i1data, i1data_sorted, a4.size(), x1.size());
    dgemm_("T", "N", x0.size(), a4.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a4.size());
  out()->add_block(odata, a4, x0);
}

void Task239::Task_local::compute() {
  const Index a4 = b(0);
  const Index x0 = b(1);
  // tensor label: I70
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma197
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x3, x2, x1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x3, x2, x1)]);
        sort_indices<1,2,3,0,0,1,2,1>(i1data, i1data_sorted, a4.size(), x3.size(), x2.size(), x1.size());
        dgemm_("T", "N", x0.size(), a4.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a4.size());
  out()->add_block(odata, a4, x0);
}

void Task240::Task_local::compute() {
  const Index a4 = b(0);
  const Index x0 = b(1);
  // tensor label: I70
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma10
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, a4, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, a4, x1)]);
        sort_indices<0,1,3,2,0,1,2,1>(i1data, i1data_sorted, x3.size(), x2.size(), a4.size(), x1.size());
        dgemm_("T", "N", x0.size(), a4.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a4.size());
  out()->add_block(odata, a4, x0);
}

void Task241::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x5, c3, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x5, c3, x4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), x5.size(), c3.size(), x4.size());
      // tensor label: I387
      std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x5, x0, x4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x5, x0, x4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, a2.size(), x5.size(), x0.size(), x4.size());
      dgemm_("T", "N", c1.size()*c3.size(), a2.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, c1.size()*c3.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), a2.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task242::Task_local::compute() {
  const Index a2 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  // tensor label: I387
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, x5, x0, x4)]);
  std::fill_n(odata.get(), out()->get_size(a2, x5, x0, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, x5, x0, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x5, x0, x4), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma126
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x5, x0, x4, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x5, x0, x4, x2, x1)]);
        sort_indices<0,4,5,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x5.size(), x0.size(), x4.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a2, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a2, x2, x1)]);
        sort_indices<0,2,3,1,0,1,-1,1>(i1data, i1data_sorted, x3.size(), a2.size(), x2.size(), x1.size());
        dgemm_("T", "N", x5.size()*x0.size()*x4.size(), a2.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x0.size()*x4.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x0.size(), x4.size(), a2.size());
  out()->add_block(odata, a2, x5, x0, x4);
}

void Task243::Task_local::compute() {
  const Index a2 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  // tensor label: I387
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, x5, x0, x4)]);
  std::fill_n(odata.get(), out()->get_size(a2, x5, x0, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, x5, x0, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x5, x0, x4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma88
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x5, x0, x4, x3, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x5, x0, x4, x3, x2)]);
        sort_indices<0,4,5,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, a2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, a2)]);
        sort_indices<2,0,1,3,0,1,-1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), a2.size());
        dgemm_("T", "N", x5.size()*x0.size()*x4.size(), a2.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x0.size()*x4.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x0.size(), x4.size(), a2.size());
  out()->add_block(odata, a2, x5, x0, x4);
}

void Task244::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c4, c1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, c4, c1, a2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), c4.size(), c1.size(), a2.size());
      // tensor label: I393
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, c4, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, c4, x0, x1)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, c3.size(), c4.size(), x0.size(), x1.size());
      dgemm_("T", "N", c1.size()*a2.size(), c3.size()*x0.size(), c4.size()*x1.size(),
             1.0, i0data_sorted, c4.size()*x1.size(), i1data_sorted, c4.size()*x1.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task245::Task_local::compute() {
  const Index c3 = b(0);
  const Index c4 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I393
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, c4, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c3, c4, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, c4, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c4, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma0
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x3, x1, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x1.size(), x2.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x3, c4, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x3, c4, x2)]);
      sort_indices<1,3,0,2,0,1,4,1>(i1data, i1data_sorted, c3.size(), x3.size(), c4.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), c3.size()*c4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c3.size(), c4.size());
  out()->add_block(odata, c3, c4, x0, x1);
}

void Task246::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c1, c4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2, c1, c4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), c4.size());
      // tensor label: I396
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, c4, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, c4, x0, x1)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, c3.size(), c4.size(), x0.size(), x1.size());
      dgemm_("T", "N", a2.size()*c1.size(), c3.size()*x0.size(), c4.size()*x1.size(),
             1.0, i0data_sorted, c4.size()*x1.size(), i1data_sorted, c4.size()*x1.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), c3.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task247::Task_local::compute() {
  const Index c3 = b(0);
  const Index c4 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I396
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, c4, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c3, c4, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, c4, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c4, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma0
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x3, x1, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x1.size(), x2.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x3, c4, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x3, c4, x2)]);
      sort_indices<1,3,0,2,0,1,-2,1>(i1data, i1data_sorted, c3.size(), x3.size(), c4.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), c3.size()*c4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c3.size(), c4.size());
  out()->add_block(odata, c3, c4, x0, x1);
}

void Task248::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c4, c3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, c4, c3, a2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), c4.size(), c3.size(), a2.size());
      // tensor label: I399
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c4, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c4, x0, x1)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), c4.size(), x0.size(), x1.size());
      dgemm_("T", "N", c3.size()*a2.size(), c1.size()*x0.size(), c4.size()*x1.size(),
             1.0, i0data_sorted, c4.size()*x1.size(), i1data_sorted, c4.size()*x1.size(),
             1.0, odata_sorted, c3.size()*a2.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task249::Task_local::compute() {
  const Index c1 = b(0);
  const Index c4 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I399
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c4, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, c4, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c4, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c4, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma0
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x3, x1, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x1.size(), x2.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3, c4, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x3, c4, x2)]);
      sort_indices<1,3,0,2,0,1,-2,1>(i1data, i1data_sorted, c1.size(), x3.size(), c4.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), c1.size()*c4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c1.size(), c4.size());
  out()->add_block(odata, c1, c4, x0, x1);
}

#endif
