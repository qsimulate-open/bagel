//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_tasks10.cc
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

#include <src/smith/mrci/MRCI_tasks10.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task450::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      // tensor label: Gamma573
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x1.size(), x0.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, x7, x6);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, x7, x6)]);
      sort_indices<2,3,0,1,0,1,1,2>(i1data, i1data_sorted, c2.size(), a1.size(), x7.size(), x6.size());
      dgemm_("T", "N", x1.size()*x0.size(), c2.size()*a1.size(), x7.size()*x6.size(),
             1.0, i0data_sorted, x7.size()*x6.size(), i1data_sorted, x7.size()*x6.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a1.size());
  out()->add_block(odata, c2, x1, x0, a1);
}

void Task451::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1, c1, a2), 0.0);
  {
    // tensor label: I108
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x1, x0, a2);
    sort_indices<2,1,0,3,1,1,1,1>(i0data, odata, c1.size(), x1.size(), x0.size(), a2.size());
  }
  out()->add_block(odata, x0, x1, c1, a2);
}

void Task452::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), a2.size());
    // tensor label: I109
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x2, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x2, x1, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), x2.size(), x1.size(), x0.size());
    dgemm_("T", "N", a2.size(), c1.size()*x1.size()*x0.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x1.size(), x0.size());
  out()->add_block(odata, c1, x1, x0, a2);
}

void Task453::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I109
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma4
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x2, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x2, x3, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, c1, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, c1, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), c1.size(), x3.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task454::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: Gamma5
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x3, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
      // tensor label: I112
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a2, c1, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a2, c1, x3)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), a2.size(), c1.size(), x3.size());
      dgemm_("T", "N", x1.size()*x0.size(), a2.size()*c1.size(), x2.size()*x3.size(),
             1.0, i0data_sorted, x2.size()*x3.size(), i1data_sorted, x2.size()*x3.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a2.size(), c1.size());
  out()->add_block(odata, c1, x1, x0, a2);
}

void Task455::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I112
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c1, x3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c1.size(), x3.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c3)]);
    sort_indices<1,0,0,1,-1,1>(i1data, i1data_sorted, x2.size(), c3.size());
    dgemm_("T", "N", a2.size()*c1.size()*x3.size(), x2.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, a2.size()*c1.size()*x3.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x3.size(), x2.size());
  out()->add_block(odata, x2, a2, c1, x3);
}

void Task456::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I112
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x3.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c3)]);
    sort_indices<1,0,0,1,2,1>(i1data, i1data_sorted, x2.size(), c3.size());
    dgemm_("T", "N", c1.size()*a2.size()*x3.size(), x2.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, c1.size()*a2.size()*x3.size());
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x3.size(), x2.size());
  out()->add_block(odata, x2, a2, c1, x3);
}

void Task457::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I112
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c3, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a2, c3, x3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c3.size(), x3.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c4, c1, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c4, c1, c3)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), c4.size(), c1.size(), c3.size());
      dgemm_("T", "N", a2.size()*x3.size(), x2.size()*c1.size(), c4.size()*c3.size(),
             1.0, i0data_sorted, c4.size()*c3.size(), i1data_sorted, c4.size()*c3.size(),
             1.0, odata_sorted, a2.size()*x3.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, a2.size(), x3.size(), x2.size(), c1.size());
  out()->add_block(odata, x2, a2, c1, x3);
}

void Task458::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I112
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c4, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c4, x3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c4.size(), x3.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c4, c1, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c4, c1, c3)]);
      sort_indices<3,1,0,2,0,1,-2,1>(i1data, i1data_sorted, x2.size(), c4.size(), c1.size(), c3.size());
      dgemm_("T", "N", a2.size()*x3.size(), x2.size()*c1.size(), c4.size()*c3.size(),
             1.0, i0data_sorted, c4.size()*c3.size(), i1data_sorted, c4.size()*c3.size(),
             1.0, odata_sorted, a2.size()*x3.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, a2.size(), x3.size(), x2.size(), c1.size());
  out()->add_block(odata, x2, a2, c1, x3);
}

void Task459::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I112
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a3, c1, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a3, c1, x3)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a3.size(), c1.size(), x3.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c4, a3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c4, a3, a2)]);
      sort_indices<1,2,0,3,0,1,-1,1>(i1data, i1data_sorted, x2.size(), c4.size(), a3.size(), a2.size());
      dgemm_("T", "N", c1.size()*x3.size(), x2.size()*a2.size(), c4.size()*a3.size(),
             1.0, i0data_sorted, c4.size()*a3.size(), i1data_sorted, c4.size()*a3.size(),
             1.0, odata_sorted, c1.size()*x3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size(), x2.size(), a2.size());
  out()->add_block(odata, x2, a2, c1, x3);
}

void Task460::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I112
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, c1, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4, c1, x3)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), c1.size(), x3.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a2, a4, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a2, a4, c3)]);
      sort_indices<3,2,0,1,0,1,2,1>(i1data, i1data_sorted, x2.size(), a2.size(), a4.size(), c3.size());
      dgemm_("T", "N", c1.size()*x3.size(), x2.size()*a2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*x3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size(), x2.size(), a2.size());
  out()->add_block(odata, x2, a2, c1, x3);
}

void Task461::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I112
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, c4, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, c4, x3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), c4.size(), x3.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c4, a3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c4, a3, a2)]);
      sort_indices<2,1,0,3,0,1,2,1>(i1data, i1data_sorted, x2.size(), c4.size(), a3.size(), a2.size());
      dgemm_("T", "N", c1.size()*x3.size(), x2.size()*a2.size(), c4.size()*a3.size(),
             1.0, i0data_sorted, c4.size()*a3.size(), i1data_sorted, c4.size()*a3.size(),
             1.0, odata_sorted, c1.size()*x3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size(), x2.size(), a2.size());
  out()->add_block(odata, x2, a2, c1, x3);
}

void Task462::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I112
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, x3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x3.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a2, a4, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a2, a4, c3)]);
      sort_indices<2,3,0,1,0,1,-1,1>(i1data, i1data_sorted, x2.size(), a2.size(), a4.size(), c3.size());
      dgemm_("T", "N", c1.size()*x3.size(), x2.size()*a2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*x3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size(), x2.size(), a2.size());
  out()->add_block(odata, x2, a2, c1, x3);
}

void Task463::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I112
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x3, x2, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x3, x2, c3)]);
      sort_indices<0,3,1,2,0,1,-1,1>(i1data, i1data_sorted, a4.size(), x3.size(), x2.size(), c3.size());
      dgemm_("T", "N", c1.size()*a2.size(), x3.size()*x2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x3.size(), x2.size());
  out()->add_block(odata, x2, a2, c1, x3);
}

void Task464::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I112
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x3, x2, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x3, x2, c3)]);
      sort_indices<3,0,1,2,0,1,2,1>(i1data, i1data_sorted, a4.size(), x3.size(), x2.size(), c3.size());
      dgemm_("T", "N", c1.size()*a2.size(), x3.size()*x2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x3.size(), x2.size());
  out()->add_block(odata, x2, a2, c1, x3);
}

void Task465::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma29
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I118
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3, a2, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x3, a2, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), x3.size(), a2.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), c1.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), a2.size());
  out()->add_block(odata, c1, x1, x0, a2);
}

void Task466::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c3, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c3, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c3.size(), x2.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), c3.size());
    dgemm_("T", "N", x3.size()*a2.size()*x2.size(), c1.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, x3.size()*a2.size()*x2.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x3.size(), a2.size(), x2.size(), c1.size());
  out()->add_block(odata, c1, x3, a2, x2);
}

void Task467::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c1, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c1.size(), x2.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, a2)]);
    sort_indices<0,1,0,1,-1,1>(i1data, i1data_sorted, a3.size(), a2.size());
    dgemm_("T", "N", x3.size()*c1.size()*x2.size(), a2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x3.size()*c1.size()*x2.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x3.size(), c1.size(), x2.size(), a2.size());
  out()->add_block(odata, c1, x3, a2, x2);
}

void Task468::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, x3, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, x3, x2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), x3.size(), x2.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c3)]);
    sort_indices<1,0,0,1,-2,1>(i1data, i1data_sorted, c1.size(), c3.size());
    dgemm_("T", "N", a2.size()*x3.size()*x2.size(), c1.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, a2.size()*x3.size()*x2.size());
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, a2.size(), x3.size(), x2.size(), c1.size());
  out()->add_block(odata, c1, x3, a2, x2);
}

void Task469::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, x3, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, x3, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), x3.size(), x2.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, a2)]);
    sort_indices<0,1,0,1,2,1>(i1data, i1data_sorted, a3.size(), a2.size());
    dgemm_("T", "N", c1.size()*x3.size()*x2.size(), a2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, c1.size()*x3.size()*x2.size());
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size(), x2.size(), a2.size());
  out()->add_block(odata, c1, x3, a2, x2);
}

void Task470::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c1, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c1.size(), a2.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x2)]);
    sort_indices<0,1,0,1,2,1>(i1data, i1data_sorted, a3.size(), x2.size());
    dgemm_("T", "N", x3.size()*c1.size()*a2.size(), x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x3.size()*c1.size()*a2.size());
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, x3.size(), c1.size(), a2.size(), x2.size());
  out()->add_block(odata, c1, x3, a2, x2);
}

void Task471::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c1, a3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c1.size(), a3.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x2)]);
    sort_indices<0,1,0,1,-1,1>(i1data, i1data_sorted, a3.size(), x2.size());
    dgemm_("T", "N", x3.size()*a2.size()*c1.size(), x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x3.size()*a2.size()*c1.size());
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x3.size(), a2.size(), c1.size(), x2.size());
  out()->add_block(odata, c1, x3, a2, x2);
}

void Task472::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c4, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c4, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c4.size(), x2.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c4, a3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c4, a3, a2)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), c4.size(), a3.size(), a2.size());
      dgemm_("T", "N", x3.size()*x2.size(), c1.size()*a2.size(), c4.size()*a3.size(),
             1.0, i0data_sorted, c4.size()*a3.size(), i1data_sorted, c4.size()*a3.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c1.size(), a2.size());
  out()->add_block(odata, c1, x3, a2, x2);
}

void Task473::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a4, c3, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c3.size(), x2.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, a4, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, a4, c3)]);
      sort_indices<2,3,0,1,0,1,-2,1>(i1data, i1data_sorted, c1.size(), a2.size(), a4.size(), c3.size());
      dgemm_("T", "N", x3.size()*x2.size(), c1.size()*a2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c1.size(), a2.size());
  out()->add_block(odata, c1, x3, a2, x2);
}

void Task474::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a3, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a3, x3, x2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a3.size(), x3.size(), x2.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c4, a3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c4, a3, a2)]);
      sort_indices<1,2,0,3,0,1,-2,1>(i1data, i1data_sorted, c1.size(), c4.size(), a3.size(), a2.size());
      dgemm_("T", "N", x3.size()*x2.size(), c1.size()*a2.size(), c4.size()*a3.size(),
             1.0, i0data_sorted, c4.size()*a3.size(), i1data_sorted, c4.size()*a3.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c1.size(), a2.size());
  out()->add_block(odata, c1, x3, a2, x2);
}

void Task475::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4, x3, x2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), x3.size(), x2.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, a4, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, a4, c3)]);
      sort_indices<3,2,0,1,0,1,4,1>(i1data, i1data_sorted, c1.size(), a2.size(), a4.size(), c3.size());
      dgemm_("T", "N", x3.size()*x2.size(), c1.size()*a2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c1.size(), a2.size());
  out()->add_block(odata, c1, x3, a2, x2);
}

void Task476::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
      // tensor label: I958
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3, x3, x2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size(), x3.size(), x2.size());
      dgemm_("T", "N", c1.size()*a2.size(), x3.size()*x2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,2,1,3,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x3.size(), x2.size());
  out()->add_block(odata, c1, x3, a2, x2);
}

void Task477::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I958
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, c3, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(a4, c3, x3, x2), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3, x3, x2);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, a4.size(), c3.size(), x3.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x3, x2, a4, c3);
    sort_indices<2,3,0,1,1,1,-2,1>(i1data, odata, x3.size(), x2.size(), a4.size(), c3.size());
  }
  out()->add_block(odata, a4, c3, x3, x2);
}

void Task478::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
      // tensor label: I961
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3, x3, x2)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size(), x3.size(), x2.size());
      dgemm_("T", "N", c1.size()*a2.size(), x3.size()*x2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,2,1,3,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x3.size(), x2.size());
  out()->add_block(odata, c1, x3, a2, x2);
}

void Task479::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I961
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, c3, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(a4, c3, x3, x2), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3, x3, x2);
    sort_indices<0,1,2,3,1,1,4,1>(i0data, odata, a4.size(), c3.size(), x3.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x3, x2, a4, c3);
    sort_indices<2,3,0,1,1,1,4,1>(i1data, odata, x3.size(), x2.size(), a4.size(), c3.size());
  }
  out()->add_block(odata, a4, c3, x3, x2);
}

void Task480::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, c4, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, c4, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), c4.size(), a2.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c4, a3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c4, a3, x2)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), c4.size(), a3.size(), x2.size());
      dgemm_("T", "N", c1.size()*a2.size(), x3.size()*x2.size(), c4.size()*a3.size(),
             1.0, i0data_sorted, c4.size()*a3.size(), i1data_sorted, c4.size()*a3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,2,1,3,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x3.size(), x2.size());
  out()->add_block(odata, c1, x3, a2, x2);
}

void Task481::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c4, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c4, a3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c4.size(), a3.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c4, a3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c4, a3, x2)]);
      sort_indices<1,2,0,3,0,1,-2,1>(i1data, i1data_sorted, x3.size(), c4.size(), a3.size(), x2.size());
      dgemm_("T", "N", c1.size()*a2.size(), x3.size()*x2.size(), c4.size()*a3.size(),
             1.0, i0data_sorted, c4.size()*a3.size(), i1data_sorted, c4.size()*a3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,2,1,3,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x3.size(), x2.size());
  out()->add_block(odata, c1, x3, a2, x2);
}

void Task482::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a4, c3, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c3.size(), a2.size());
      // tensor label: I1006
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x2, c1, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x2, c1, c3)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, a4.size(), x2.size(), c1.size(), c3.size());
      dgemm_("T", "N", x3.size()*a2.size(), x2.size()*c1.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, x3.size()*a2.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x3.size(), a2.size(), x2.size(), c1.size());
  out()->add_block(odata, c1, x3, a2, x2);
}

void Task483::Task_local::compute() {
  const Index a4 = b(0);
  const Index x2 = b(1);
  const Index c1 = b(2);
  const Index c3 = b(3);
  // tensor label: I1006
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, x2, c1, c3)]);
  std::fill_n(odata.get(), out()->get_size(a4, x2, c1, c3), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x2, c1, c3);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, a4.size(), x2.size(), c1.size(), c3.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, x2, a4, c3);
    sort_indices<2,1,0,3,1,1,1,1>(i1data, odata, c1.size(), x2.size(), a4.size(), c3.size());
  }
  out()->add_block(odata, a4, x2, c1, c3);
}

void Task484::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c3, a4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c3.size(), a4.size());
      // tensor label: I1009
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x2, c1, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x2, c1, c3)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a4.size(), x2.size(), c1.size(), c3.size());
      dgemm_("T", "N", x3.size()*a2.size(), x2.size()*c1.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, x3.size()*a2.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x3.size(), a2.size(), x2.size(), c1.size());
  out()->add_block(odata, c1, x3, a2, x2);
}

void Task485::Task_local::compute() {
  const Index a4 = b(0);
  const Index x2 = b(1);
  const Index c1 = b(2);
  const Index c3 = b(3);
  // tensor label: I1009
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, x2, c1, c3)]);
  std::fill_n(odata.get(), out()->get_size(a4, x2, c1, c3), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x2, c1, c3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, a4.size(), x2.size(), c1.size(), c3.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, x2, a4, c3);
    sort_indices<2,1,0,3,1,1,-2,1>(i1data, odata, c1.size(), x2.size(), a4.size(), c3.size());
  }
  out()->add_block(odata, a4, x2, c1, c3);
}

void Task486::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c1, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a4, c1, a3)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c1.size(), a3.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x2, a3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x2, a3, a2)]);
      sort_indices<0,2,1,3,0,1,2,1>(i1data, i1data_sorted, a4.size(), x2.size(), a3.size(), a2.size());
      dgemm_("T", "N", x3.size()*c1.size(), x2.size()*a2.size(), a4.size()*a3.size(),
             1.0, i0data_sorted, a4.size()*a3.size(), i1data_sorted, a4.size()*a3.size(),
             1.0, odata_sorted, x3.size()*c1.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x3.size(), c1.size(), x2.size(), a2.size());
  out()->add_block(odata, c1, x3, a2, x2);
}

void Task487::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c1, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c1, a4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c1.size(), a4.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x2, a3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x2, a3, a2)]);
      sort_indices<2,0,1,3,0,1,-1,1>(i1data, i1data_sorted, a4.size(), x2.size(), a3.size(), a2.size());
      dgemm_("T", "N", x3.size()*c1.size(), x2.size()*a2.size(), a4.size()*a3.size(),
             1.0, i0data_sorted, a4.size()*a3.size(), i1data_sorted, a4.size()*a3.size(),
             1.0, odata_sorted, x3.size()*c1.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x3.size(), c1.size(), x2.size(), a2.size());
  out()->add_block(odata, c1, x3, a2, x2);
}

void Task488::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), x2.size());
    // tensor label: I130
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x2, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x2, x1, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), x2.size(), x1.size(), x0.size());
    dgemm_("T", "N", c1.size(), a2.size()*x1.size()*x0.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, c1.size());
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), x0.size());
  out()->add_block(odata, c1, x1, x0, a2);
}

void Task489::Task_local::compute() {
  const Index a2 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I130
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a2, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma252
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2, x4, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x2, x4, x3, x1, x0)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a2, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a2, x4, x3)]);
        sort_indices<0,2,3,1,0,1,-1,1>(i1data, i1data_sorted, x5.size(), a2.size(), x4.size(), x3.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), a2.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), a2.size());
  out()->add_block(odata, a2, x2, x1, x0);
}

void Task490::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  // tensor label: Gamma32
  std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
  // tensor label: I133
  std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size());
  dgemm_("T", "N", x1.size()*x0.size(), c1.size()*a2.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, x1.size()*x0.size());
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), a2.size());
  out()->add_block(odata, c1, x1, x0, a2);
}

void Task491::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: I133
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
      // tensor label: h1
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3)]);
      sort_indices<0,1,0,1,-4,1>(i1data, i1data_sorted, a4.size(), c3.size());
      dgemm_("T", "N", c1.size()*a2.size(), 1, a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size());
  out()->add_block(odata, c1, a2);
}

void Task492::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: I133
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
      // tensor label: h1
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3)]);
      sort_indices<1,0,0,1,8,1>(i1data, i1data_sorted, a4.size(), c3.size());
      dgemm_("T", "N", c1.size()*a2.size(), 1, a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size());
  out()->add_block(odata, c1, a2);
}

void Task493::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: I133
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      for (auto& c5 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, c5, a2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4, c5, a2)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), c5.size(), a2.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c5, a4, c3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c5, a4, c3)]);
        sort_indices<3,2,1,0,0,1,-8,1>(i1data, i1data_sorted, c1.size(), c5.size(), a4.size(), c3.size());
        dgemm_("T", "N", a2.size(), c1.size(), c5.size()*a4.size()*c3.size(),
               1.0, i0data_sorted, c5.size()*a4.size()*c3.size(), i1data_sorted, c5.size()*a4.size()*c3.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size());
  out()->add_block(odata, c1, a2);
}

void Task494::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: I133
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& c5 : *range_[0]) {
      for (auto& a4 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c5, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c5, a4)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c5.size(), a4.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c5, a4, c3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c5, a4, c3)]);
        sort_indices<3,1,2,0,0,1,4,1>(i1data, i1data_sorted, c1.size(), c5.size(), a4.size(), c3.size());
        dgemm_("T", "N", a2.size(), c1.size(), c5.size()*a4.size()*c3.size(),
               1.0, i0data_sorted, c5.size()*a4.size()*c3.size(), i1data_sorted, c5.size()*a4.size()*c3.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size());
  out()->add_block(odata, c1, a2);
}

void Task495::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: I133
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2), 0.0);
  for (auto& a5 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      for (auto& a4 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a5, c3, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a5, c3, a4)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), a5.size(), c3.size(), a4.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(a5, a2, a4, c3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, a2, a4, c3)]);
        sort_indices<0,3,2,1,0,1,8,1>(i1data, i1data_sorted, a5.size(), a2.size(), a4.size(), c3.size());
        dgemm_("T", "N", c1.size(), a2.size(), a5.size()*a4.size()*c3.size(),
               1.0, i0data_sorted, a5.size()*a4.size()*c3.size(), i1data_sorted, a5.size()*a4.size()*c3.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size());
  out()->add_block(odata, c1, a2);
}

void Task496::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: I133
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      for (auto& a5 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a5);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a5)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a5.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(a5, a2, a4, c3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, a2, a4, c3)]);
        sort_indices<2,3,0,1,0,1,-4,1>(i1data, i1data_sorted, a5.size(), a2.size(), a4.size(), c3.size());
        dgemm_("T", "N", c1.size(), a2.size(), a5.size()*a4.size()*c3.size(),
               1.0, i0data_sorted, a5.size()*a4.size()*c3.size(), i1data_sorted, a5.size()*a4.size()*c3.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size());
  out()->add_block(odata, c1, a2);
}

void Task497::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, x2, c3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, x2, c3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), x2.size(), c3.size());
        // tensor label: I840
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c3, x3, x2, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c3, x3, x2, x1, x0)]);
        sort_indices<2,3,1,0,4,5,0,1,1,1>(i1data, i1data_sorted, c1.size(), c3.size(), x3.size(), x2.size(), x1.size(), x0.size());
        dgemm_("T", "N", a2.size(), c1.size()*x1.size()*x0.size(), c3.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, c3.size()*x3.size()*x2.size(), i1data_sorted, c3.size()*x3.size()*x2.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x1.size(), x0.size());
  out()->add_block(odata, c1, x1, x0, a2);
}

void Task498::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I840
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x3, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x3, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x3, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: Gamma107
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x5, x2, x4, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x5, x2, x4, x1, x0)]);
      sort_indices<1,3,0,2,4,5,0,1,1,1>(i0data, i0data_sorted, x3.size(), x5.size(), x2.size(), x4.size(), x1.size(), x0.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x5, c3, x4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x5, c3, x4)]);
      sort_indices<1,3,0,2,0,1,-2,1>(i1data, i1data_sorted, c1.size(), x5.size(), c3.size(), x4.size());
      dgemm_("T", "N", x3.size()*x2.size()*x1.size()*x0.size(), c1.size()*c3.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, x3.size()*x2.size()*x1.size()*x0.size());
    }
  }
  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size(), c1.size(), c3.size());
  out()->add_block(odata, c1, c3, x3, x2, x1, x0);
}

void Task499::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x4, a2, x3, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x4, a2, x3, x2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x4.size(), a2.size(), x3.size(), x2.size());
        // tensor label: I843
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x4, x3, x2, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x4, x3, x2, x1, x0)]);
        sort_indices<1,2,3,0,4,5,0,1,1,1>(i1data, i1data_sorted, c1.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());
        dgemm_("T", "N", a2.size(), c1.size()*x1.size()*x0.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x1.size(), x0.size());
  out()->add_block(odata, c1, x1, x0, a2);
}

#endif
