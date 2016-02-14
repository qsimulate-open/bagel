//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_tasks18.cc
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

#include <src/smith/mrci/MRCI_tasks18.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task850::Task_local::compute() {
  const Index a4 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1508
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, a1, x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, a1, x3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a1, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a1, x3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma503
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x2, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x2, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x2.size(), x0.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a1, a4, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a1, a4, x1)]);
      sort_indices<3,0,1,2,0,1,-1,1>(i1data, i1data_sorted, x2.size(), a1.size(), a4.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), a1.size()*a4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<3,2,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a1.size(), a4.size());
  out()->add_block(odata, a4, a1, x3, x0);
}

void Task851::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a4.size());
      // tensor label: I1511
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a1, x3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a1, x3, x0)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), a1.size(), x3.size(), x0.size());
      dgemm_("T", "N", a3.size()*c2.size(), a1.size()*x0.size(), a4.size()*x3.size(),
             1.0, i0data_sorted, a4.size()*x3.size(), i1data_sorted, a4.size()*x3.size(),
             1.0, odata_sorted, a3.size()*c2.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task852::Task_local::compute() {
  const Index a4 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1511
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, a1, x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, a1, x3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a1, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a1, x3, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma51
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I1512
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a1, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a1, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a1.size(), x2.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), a4.size()*a1.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a4.size(), a1.size());
  out()->add_block(odata, a4, a1, x3, x0);
}

void Task853::Task_local::compute() {
  const Index a4 = b(0);
  const Index a1 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1512
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, a1, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a4, a1, x2, x1), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, a4.size(), a1.size(), x2.size(), x1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x2, x1, a4, a1);
    sort_indices<2,3,0,1,1,1,-1,2>(i1data, odata, x2.size(), x1.size(), a4.size(), a1.size());
  }
  out()->add_block(odata, a4, a1, x2, x1);
}

void Task854::Task_local::compute() {
  const Index a4 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1511
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, a1, x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, a1, x3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a1, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a1, x3, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma29
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x2, x1, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x2, x1, a1)]);
      sort_indices<1,2,0,3,0,1,-1,2>(i1data, i1data_sorted, a4.size(), x2.size(), x1.size(), a1.size());
      dgemm_("T", "N", x3.size()*x0.size(), a4.size()*a1.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a4.size(), a1.size());
  out()->add_block(odata, a4, a1, x3, x0);
}

void Task855::Task_local::compute() {
  const Index a4 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1511
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, a1, x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, a1, x3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a1, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a1, x3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma503
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x2, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x2, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x2.size(), x0.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a1, a4, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a1, a4, x1)]);
      sort_indices<3,0,1,2,0,1,1,2>(i1data, i1data_sorted, x2.size(), a1.size(), a4.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), a1.size()*a4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<3,2,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a1.size(), a4.size());
  out()->add_block(odata, a4, a1, x3, x0);
}

void Task856::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a4, c2, a1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c2.size(), a1.size());
      // tensor label: I1514
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a3, x3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a3, x3, x0)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), a3.size(), x3.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), a3.size()*x0.size(), a4.size()*x3.size(),
             1.0, i0data_sorted, a4.size()*x3.size(), i1data_sorted, a4.size()*x3.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task857::Task_local::compute() {
  const Index a4 = b(0);
  const Index a3 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1514
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, a3, x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, a3, x3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a3, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a3, x3, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma51
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I1515
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a3, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a3, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a3.size(), x2.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), a4.size()*a3.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a4.size(), a3.size());
  out()->add_block(odata, a4, a3, x3, x0);
}

void Task858::Task_local::compute() {
  const Index a4 = b(0);
  const Index a3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1515
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, a3, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a4, a3, x2, x1), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a3, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, a4.size(), a3.size(), x2.size(), x1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x2, x1, a4, a3);
    sort_indices<2,3,0,1,1,1,-1,2>(i1data, odata, x2.size(), x1.size(), a4.size(), a3.size());
  }
  out()->add_block(odata, a4, a3, x2, x1);
}

void Task859::Task_local::compute() {
  const Index a4 = b(0);
  const Index a3 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1514
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, a3, x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, a3, x3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a3, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a3, x3, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma29
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x2, x1, a3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x2, x1, a3)]);
      sort_indices<1,2,0,3,0,1,-1,2>(i1data, i1data_sorted, a4.size(), x2.size(), x1.size(), a3.size());
      dgemm_("T", "N", x3.size()*x0.size(), a4.size()*a3.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a4.size(), a3.size());
  out()->add_block(odata, a4, a3, x3, x0);
}

void Task860::Task_local::compute() {
  const Index a4 = b(0);
  const Index a3 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1514
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, a3, x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, a3, x3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a3, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a3, x3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma503
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x2, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x2, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x2.size(), x0.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a3, a4, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a3, a4, x1)]);
      sort_indices<3,0,1,2,0,1,1,2>(i1data, i1data_sorted, x2.size(), a3.size(), a4.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), a3.size()*a4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<3,2,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a3.size(), a4.size());
  out()->add_block(odata, a4, a3, x3, x0);
}

void Task861::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a4.size());
      // tensor label: I1517
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a3, x3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a3, x3, x0)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), a3.size(), x3.size(), x0.size());
      dgemm_("T", "N", a1.size()*c2.size(), a3.size()*x0.size(), a4.size()*x3.size(),
             1.0, i0data_sorted, a4.size()*x3.size(), i1data_sorted, a4.size()*x3.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<2,1,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task862::Task_local::compute() {
  const Index a4 = b(0);
  const Index a3 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1517
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, a3, x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, a3, x3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a3, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a3, x3, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma51
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I1518
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a3, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a3, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a3.size(), x2.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), a4.size()*a3.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a4.size(), a3.size());
  out()->add_block(odata, a4, a3, x3, x0);
}

void Task863::Task_local::compute() {
  const Index a4 = b(0);
  const Index a3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1518
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, a3, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a4, a3, x2, x1), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a3, x2, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, a4.size(), a3.size(), x2.size(), x1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x2, a3, a4, x1);
    sort_indices<2,1,0,3,1,1,-1,2>(i1data, odata, x2.size(), a3.size(), a4.size(), x1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i2data = in(0)->get_block(x2, x1, a4, a3);
    sort_indices<2,3,0,1,1,1,1,1>(i2data, odata, x2.size(), x1.size(), a4.size(), a3.size());
  }
  out()->add_block(odata, a4, a3, x2, x1);
}

void Task864::Task_local::compute() {
  const Index a4 = b(0);
  const Index a3 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1517
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, a3, x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, a3, x3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a3, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a3, x3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma27
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x1, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x2, x1, a3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x2, x1, a3)]);
      sort_indices<2,1,0,3,0,1,1,2>(i1data, i1data_sorted, a4.size(), x2.size(), x1.size(), a3.size());
      dgemm_("T", "N", x3.size()*x0.size(), a4.size()*a3.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a4.size(), a3.size());
  out()->add_block(odata, a4, a3, x3, x0);
}

void Task865::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a3.size());
      // tensor label: I1604
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x5, x0, x4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x5, x0, x4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, c2.size(), x5.size(), x0.size(), x4.size());
      dgemm_("T", "N", a1.size()*a3.size(), c2.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, a1.size()*a3.size());
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), c2.size(), x0.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task866::Task_local::compute() {
  const Index c2 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  // tensor label: I1604
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x5, x0, x4)]);
  std::fill_n(odata.get(), out()->get_size(c2, x5, x0, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x5, x0, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x5, x0, x4), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma50
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<3,4,5,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x3, x2, x1)]);
        sort_indices<1,2,3,0,0,1,-1,1>(i1data, i1data_sorted, c2.size(), x3.size(), x2.size(), x1.size());
        dgemm_("T", "N", x5.size()*x0.size()*x4.size(), c2.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x0.size()*x4.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x0.size(), x4.size(), c2.size());
  out()->add_block(odata, c2, x5, x0, x4);
}

void Task867::Task_local::compute() {
  const Index c2 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  // tensor label: I1604
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x5, x0, x4)]);
  std::fill_n(odata.get(), out()->get_size(c2, x5, x0, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x5, x0, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x5, x0, x4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma526
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x1, x3, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x1, x3, x2)]);
        sort_indices<3,4,5,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x1.size(), x3.size(), x2.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, c2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, c2, x1)]);
        sort_indices<3,0,1,2,0,1,-1,1>(i1data, i1data_sorted, x3.size(), x2.size(), c2.size(), x1.size());
        dgemm_("T", "N", x5.size()*x0.size()*x4.size(), c2.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x0.size()*x4.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x0.size(), x4.size(), c2.size());
  out()->add_block(odata, c2, x5, x0, x4);
}

void Task868::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x1, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, x1, c2, a1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, a4.size(), x1.size(), c2.size(), a1.size());
      // tensor label: I1610
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, a4, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, a4, x0, x1)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, a3.size(), a4.size(), x0.size(), x1.size());
      dgemm_("T", "N", c2.size()*a1.size(), a3.size()*x0.size(), a4.size()*x1.size(),
             1.0, i0data_sorted, a4.size()*x1.size(), i1data_sorted, a4.size()*x1.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task869::Task_local::compute() {
  const Index a3 = b(0);
  const Index a4 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1610
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, a4, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a3, a4, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, a4, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, a4, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma51
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a3, x2, a4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a3, x2, a4)]);
      sort_indices<0,2,1,3,0,1,-2,1>(i1data, i1data_sorted, x3.size(), a3.size(), x2.size(), a4.size());
      dgemm_("T", "N", x0.size()*x1.size(), a3.size()*a4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a3.size(), a4.size());
  out()->add_block(odata, a3, a4, x0, x1);
}

void Task870::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x1, c2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, x1, c2, a3)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, a4.size(), x1.size(), c2.size(), a3.size());
      // tensor label: I1613
      std::unique_ptr<double[]> i1data = in(1)->get_block(a1, a4, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, a4, x0, x1)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, a1.size(), a4.size(), x0.size(), x1.size());
      dgemm_("T", "N", c2.size()*a3.size(), a1.size()*x0.size(), a4.size()*x1.size(),
             1.0, i0data_sorted, a4.size()*x1.size(), i1data_sorted, a4.size()*x1.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), a1.size(), x0.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task871::Task_local::compute() {
  const Index a1 = b(0);
  const Index a4 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1613
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, a4, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, a4, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, a4, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, a4, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma51
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, x2, a4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, x2, a4)]);
      sort_indices<0,2,1,3,0,1,4,1>(i1data, i1data_sorted, x3.size(), a1.size(), x2.size(), a4.size());
      dgemm_("T", "N", x0.size()*x1.size(), a1.size()*a4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), a4.size());
  out()->add_block(odata, a1, a4, x0, x1);
}

void Task872::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x1, a4, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x1, a4, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), x1.size(), a4.size(), a1.size());
      // tensor label: I1616
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, a4, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, a4, x1, x0)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), a4.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), a3.size()*x0.size(), a4.size()*x1.size(),
             1.0, i0data_sorted, a4.size()*x1.size(), i1data_sorted, a4.size()*x1.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task873::Task_local::compute() {
  const Index a3 = b(0);
  const Index a4 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1616
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, a4, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a3, a4, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, a4, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, a4, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma503
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x2, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x2, x0)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x2.size(), x0.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a3, x2, a4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a3, x2, a4)]);
      sort_indices<0,2,1,3,0,1,2,1>(i1data, i1data_sorted, x3.size(), a3.size(), x2.size(), a4.size());
      dgemm_("T", "N", x1.size()*x0.size(), a3.size()*a4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a3.size(), a4.size());
  out()->add_block(odata, a3, a4, x1, x0);
}

void Task874::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x1, a4, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x1, a4, a3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), x1.size(), a4.size(), a3.size());
      // tensor label: I1619
      std::unique_ptr<double[]> i1data = in(1)->get_block(a1, a4, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, a4, x0, x1)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, a1.size(), a4.size(), x0.size(), x1.size());
      dgemm_("T", "N", c2.size()*a3.size(), a1.size()*x0.size(), a4.size()*x1.size(),
             1.0, i0data_sorted, a4.size()*x1.size(), i1data_sorted, a4.size()*x1.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), a1.size(), x0.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task875::Task_local::compute() {
  const Index a1 = b(0);
  const Index a4 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1619
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, a4, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, a4, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, a4, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, a4, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma51
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, x2, a4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, x2, a4)]);
      sort_indices<0,2,1,3,0,1,-2,1>(i1data, i1data_sorted, x3.size(), a1.size(), x2.size(), a4.size());
      dgemm_("T", "N", x0.size()*x1.size(), a1.size()*a4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), a4.size());
  out()->add_block(odata, a1, a4, x0, x1);
}

void Task876::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: Gamma562
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size());
    // tensor label: I1701
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a3, c2, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a3.size(), c2.size(), a1.size());
    dgemm_("T", "N", x0.size(), a3.size()*c2.size()*a1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,2,0,3,1,1,1,1>(odata_sorted, odata, x0.size(), a3.size(), c2.size(), a1.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task877::Task_local::compute() {
  const Index x3 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I1701
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a3, c2, a1)]);
  std::fill_n(odata.get(), out()->get_size(x3, a3, c2, a1), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), a3.size(), c2.size(), a1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x3, a1, c2, a3);
    sort_indices<0,3,2,1,1,1,-2,1>(i1data, odata, x3.size(), a1.size(), c2.size(), a3.size());
  }
  out()->add_block(odata, x3, a3, c2, a1);
}

void Task878::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    // tensor label: Gamma564
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size());
    // tensor label: I1705
    std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a3, c2, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), a3.size(), c2.size(), a1.size());
    dgemm_("T", "N", x0.size(), a3.size()*c2.size()*a1.size(), x5.size(),
           1.0, i0data_sorted, x5.size(), i1data_sorted, x5.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,2,0,3,1,1,1,1>(odata_sorted, odata, x0.size(), a3.size(), c2.size(), a1.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task879::Task_local::compute() {
  const Index x5 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I1705
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, a3, c2, a1)]);
  std::fill_n(odata.get(), out()->get_size(x5, a3, c2, a1), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, c2, a1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x5.size(), a3.size(), c2.size(), a1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x5, a1, c2, a3);
    sort_indices<0,3,2,1,1,1,-1,1>(i1data, odata, x5.size(), a1.size(), c2.size(), a3.size());
  }
  out()->add_block(odata, x5, a3, c2, a1);
}

void Task880::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, a2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(x1, a2, x0, a1), 0.0);
  {
    // tensor label: I239
    std::unique_ptr<double[]> i0data = in(0)->get_block(a1, x0, x1, a2);
    sort_indices<2,3,1,0,1,1,1,1>(i0data, odata, a1.size(), x0.size(), x1.size(), a2.size());
  }
  {
    // tensor label: I239
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, x1, x0, a1);
    sort_indices<1,0,2,3,1,1,1,1>(i0data, odata, a2.size(), x1.size(), x0.size(), a1.size());
  }
  out()->add_block(odata, x1, a2, x0, a1);
}

void Task881::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  // tensor label: I239
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), a2.size());
    // tensor label: I240
    std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x0, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x0, x2, x1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a1.size(), x0.size(), x2.size(), x1.size());
    dgemm_("T", "N", a2.size(), a1.size()*x0.size()*x1.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), a1.size(), x0.size(), x1.size());
  out()->add_block(odata, a1, x0, x1, a2);
}

void Task882::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I240
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma50
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: t2
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
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task883::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  // tensor label: I239
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma51
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I243
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x3, a1, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x3, a1, a2)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), x3.size(), a1.size(), a2.size());
      dgemm_("T", "N", x0.size()*x1.size(), a1.size()*a2.size(), x2.size()*x3.size(),
             1.0, i0data_sorted, x2.size()*x3.size(), i1data_sorted, x2.size()*x3.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), a2.size());
  out()->add_block(odata, a1, x0, x1, a2);
}

void Task884::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I243
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, x3, a1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x2, x3, a1, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x3, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, a1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c3, a2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c3.size(), a2.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c3)]);
    sort_indices<1,0,0,1,-1,1>(i1data, i1data_sorted, x2.size(), c3.size());
    dgemm_("T", "N", x3.size()*a1.size()*a2.size(), x2.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, x3.size()*a1.size()*a2.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size(), a2.size(), x2.size());
  out()->add_block(odata, x2, x3, a1, a2);
}

void Task885::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I243
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, x3, a1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x2, x3, a1, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x3, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, a1, a2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, a2)]);
    sort_indices<0,1,0,1,2,1>(i1data, i1data_sorted, a3.size(), a2.size());
    dgemm_("T", "N", x3.size()*a1.size()*x2.size(), a2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x3.size()*a1.size()*x2.size());
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size(), x2.size(), a2.size());
  out()->add_block(odata, x2, x3, a1, a2);
}

void Task886::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I243
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, x3, a1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x2, x3, a1, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x3, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, a1, a2), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c4, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c4, a3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c4.size(), a3.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c4, a3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c4, a3, a2)]);
      sort_indices<1,2,0,3,0,1,-1,1>(i1data, i1data_sorted, x2.size(), c4.size(), a3.size(), a2.size());
      dgemm_("T", "N", x3.size()*a1.size(), x2.size()*a2.size(), c4.size()*a3.size(),
             1.0, i0data_sorted, c4.size()*a3.size(), i1data_sorted, c4.size()*a3.size(),
             1.0, odata_sorted, x3.size()*a1.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size(), x2.size(), a2.size());
  out()->add_block(odata, x2, x3, a1, a2);
}

void Task887::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I243
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, x3, a1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x2, x3, a1, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x3, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, a1, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c3, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a4, c3, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c3.size(), a1.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a2, a4, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a2, a4, c3)]);
      sort_indices<2,3,0,1,0,1,-1,1>(i1data, i1data_sorted, x2.size(), a2.size(), a4.size(), c3.size());
      dgemm_("T", "N", x3.size()*a1.size(), x2.size()*a2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, x3.size()*a1.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size(), x2.size(), a2.size());
  out()->add_block(odata, x2, x3, a1, a2);
}

void Task888::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I243
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, x3, a1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x2, x3, a1, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x3, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, a1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c3, a4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c3.size(), a4.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a2, a4, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a2, a4, c3)]);
      sort_indices<3,2,0,1,0,1,2,1>(i1data, i1data_sorted, x2.size(), a2.size(), a4.size(), c3.size());
      dgemm_("T", "N", x3.size()*a1.size(), x2.size()*a2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, x3.size()*a1.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size(), x2.size(), a2.size());
  out()->add_block(odata, x2, x3, a1, a2);
}

void Task889::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  // tensor label: I239
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& c3 : *range_[0]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c3, x4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c3, x4)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c3.size(), x4.size());
        // tensor label: I1622
        std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c3, x5, x0, x4, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c3, x5, x0, x4, x1)]);
        sort_indices<2,1,4,0,3,5,0,1,1,1>(i1data, i1data_sorted, a2.size(), c3.size(), x5.size(), x0.size(), x4.size(), x1.size());
        dgemm_("T", "N", a1.size(), a2.size()*x0.size()*x1.size(), c3.size()*x5.size()*x4.size(),
               1.0, i0data_sorted, c3.size()*x5.size()*x4.size(), i1data_sorted, c3.size()*x5.size()*x4.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
  out()->add_block(odata, a1, x0, x1, a2);
}

void Task890::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x5 = b(2);
  const Index x0 = b(3);
  const Index x4 = b(4);
  const Index x1 = b(5);
  // tensor label: I1622
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c3, x5, x0, x4, x1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c3, x5, x0, x4, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c3, x5, x0, x4, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x5, x0, x4, x1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: Gamma531
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x2, x4, x3, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x2, x4, x3, x1)]);
      sort_indices<2,4,0,1,3,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x2.size(), x4.size(), x3.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a2, x2, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a2, x2, c3)]);
      sort_indices<2,0,1,3,0,1,-1,1>(i1data, i1data_sorted, x3.size(), a2.size(), x2.size(), c3.size());
      dgemm_("T", "N", x5.size()*x0.size()*x4.size()*x1.size(), a2.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x5.size()*x0.size()*x4.size()*x1.size());
    }
  }
  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, x5.size(), x0.size(), x4.size(), x1.size(), a2.size(), c3.size());
  out()->add_block(odata, a2, c3, x5, x0, x4, x1);
}

void Task891::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  // tensor label: I239
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& x5 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a1, x5, x4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a1, x5, x4)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a1.size(), x5.size(), x4.size());
        // tensor label: I1625
        std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c3, x5, x4, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c3, x5, x4, x1, x0)]);
        sort_indices<1,2,3,0,4,5,0,1,1,1>(i1data, i1data_sorted, a2.size(), c3.size(), x5.size(), x4.size(), x1.size(), x0.size());
        dgemm_("T", "N", a1.size(), a2.size()*x1.size()*x0.size(), c3.size()*x5.size()*x4.size(),
               1.0, i0data_sorted, c3.size()*x5.size()*x4.size(), i1data_sorted, c3.size()*x5.size()*x4.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x1.size(), x0.size());
  out()->add_block(odata, a1, x0, x1, a2);
}

void Task892::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I1625
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c3, x5, x4, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a2, c3, x5, x4, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c3, x5, x4, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x5, x4, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma532
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x1, x2, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x1, x2, x0)]);
      sort_indices<2,4,0,1,3,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x1.size(), x2.size(), x0.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a2, x2, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a2, x2, c3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a2.size(), x2.size(), c3.size());
      dgemm_("T", "N", x5.size()*x4.size()*x1.size()*x0.size(), a2.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x5.size()*x4.size()*x1.size()*x0.size());
    }
  }
  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x1.size(), x0.size(), a2.size(), c3.size());
  out()->add_block(odata, a2, c3, x5, x4, x1, x0);
}

void Task893::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  // tensor label: I239
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a1, x6, x5);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a1, x6, x5)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x7.size(), a1.size(), x6.size(), x5.size());
        // tensor label: I1628
        std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x7, x0, x6, x5, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x7, x0, x6, x5, x1)]);
        sort_indices<1,3,4,0,2,5,0,1,1,1>(i1data, i1data_sorted, a2.size(), x7.size(), x0.size(), x6.size(), x5.size(), x1.size());
        dgemm_("T", "N", a1.size(), a2.size()*x0.size()*x1.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
  out()->add_block(odata, a1, x0, x1, a2);
}

void Task894::Task_local::compute() {
  const Index a2 = b(0);
  const Index x7 = b(1);
  const Index x0 = b(2);
  const Index x6 = b(3);
  const Index x5 = b(4);
  const Index x1 = b(5);
  // tensor label: I1628
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, x7, x0, x6, x5, x1)]);
  std::fill_n(odata.get(), out()->get_size(a2, x7, x0, x6, x5, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, x7, x0, x6, x5, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x7, x0, x6, x5, x1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma533
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x5, x4, x1, x3, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x0, x6, x5, x4, x1, x3, x2)]);
        sort_indices<4,6,7,0,1,2,3,5,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x6.size(), x5.size(), x4.size(), x1.size(), x3.size(), x2.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, a2, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, a2, x3, x2)]);
        sort_indices<0,2,3,1,0,1,1,2>(i1data, i1data_sorted, x4.size(), a2.size(), x3.size(), x2.size());
        dgemm_("T", "N", x7.size()*x0.size()*x6.size()*x5.size()*x1.size(), a2.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x7.size()*x0.size()*x6.size()*x5.size()*x1.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x7.size(), x0.size(), x6.size(), x5.size(), x1.size(), a2.size());
  out()->add_block(odata, a2, x7, x0, x6, x5, x1);
}

void Task895::Task_local::compute() {
  const Index a2 = b(0);
  const Index x7 = b(1);
  const Index x0 = b(2);
  const Index x6 = b(3);
  const Index x5 = b(4);
  const Index x1 = b(5);
  // tensor label: I1628
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, x7, x0, x6, x5, x1)]);
  std::fill_n(odata.get(), out()->get_size(a2, x7, x0, x6, x5, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, x7, x0, x6, x5, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x7, x0, x6, x5, x1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma349
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x5, x4, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x0, x6, x5, x4, x3, x2, x1)]);
        sort_indices<4,5,6,0,1,2,3,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x6.size(), x5.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, x2, a2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, x2, a2)]);
        sort_indices<0,1,2,3,0,1,1,2>(i1data, i1data_sorted, x4.size(), x3.size(), x2.size(), a2.size());
        dgemm_("T", "N", x7.size()*x0.size()*x6.size()*x5.size()*x1.size(), a2.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x7.size()*x0.size()*x6.size()*x5.size()*x1.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x7.size(), x0.size(), x6.size(), x5.size(), x1.size(), a2.size());
  out()->add_block(odata, a2, x7, x0, x6, x5, x1);
}

void Task896::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  // tensor label: I239
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a1, a3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a1, a3, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x2.size(), a1.size(), a3.size(), a2.size());
      // tensor label: I1634
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x1, x2, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x1, x2, x0)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), x1.size(), x2.size(), x0.size());
      dgemm_("T", "N", a1.size()*a2.size(), x1.size()*x0.size(), a3.size()*x2.size(),
             1.0, i0data_sorted, a3.size()*x2.size(), i1data_sorted, a3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x1.size(), x0.size());
  out()->add_block(odata, a1, x0, x1, a2);
}

void Task897::Task_local::compute() {
  const Index a3 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index x0 = b(3);
  // tensor label: I1634
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, x1, x2, x0)]);
  std::fill_n(odata.get(), out()->get_size(a3, x1, x2, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x1, x2, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x1, x2, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma471
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x1, x4, x3, x2, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x1, x4, x3, x2, x0)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x1.size(), x4.size(), x3.size(), x2.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a3, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a3, x4, x3)]);
        sort_indices<0,2,3,1,0,1,-1,1>(i1data, i1data_sorted, x5.size(), a3.size(), x4.size(), x3.size());
        dgemm_("T", "N", x1.size()*x2.size()*x0.size(), a3.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x1.size()*x2.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x2.size(), x0.size(), a3.size());
  out()->add_block(odata, a3, x1, x2, x0);
}

void Task898::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  // tensor label: I239
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c3, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c3.size(), a2.size());
      // tensor label: I1640
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x5, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x5, x0, x1)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), x5.size(), x0.size(), x1.size());
      dgemm_("T", "N", a1.size()*a2.size(), x0.size()*x1.size(), c3.size()*x5.size(),
             1.0, i0data_sorted, c3.size()*x5.size(), i1data_sorted, c3.size()*x5.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
  out()->add_block(odata, a1, x0, x1, a2);
}

void Task899::Task_local::compute() {
  const Index c3 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1640
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, x5, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c3, x5, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x5, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x5, x0, x1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma526
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x1, x3, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x1, x3, x2)]);
        sort_indices<2,4,5,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x1.size(), x3.size(), x2.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, c3, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, c3, x3, x2)]);
        sort_indices<0,2,3,1,0,1,-1,2>(i1data, i1data_sorted, x4.size(), c3.size(), x3.size(), x2.size());
        dgemm_("T", "N", x5.size()*x0.size()*x1.size(), c3.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x5.size()*x0.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x0.size(), x1.size(), c3.size());
  out()->add_block(odata, c3, x5, x0, x1);
}

#endif
